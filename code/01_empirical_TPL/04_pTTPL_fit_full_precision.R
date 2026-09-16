# ==============================================================================
# NWS pTTPL REVISION WORKFLOW — STEP 1
# Re-fit pTTPL from MPD-weighted annual abundance matrices with full precision
#
# INPUT, for example:
#   ./predata/Archaea/UPARSE/N01_MPD_weighted_otu_All_abun.csv
#   ./predata/Bacteria/UPARSE/N01_MPD_weighted_otu_All_abun.csv
#   ./predata/Fungi/UPARSE/N01_MPD_weighted_otu_All_abun.csv
#   ./predata/Protists/UPARSE/N01_MPD_weighted_otu_All_abun.csv
#   ./predata/Nematoda/UPARSE/N01_MPD_weighted_otu_All_abun.csv
#
# DESIGN:
#   5 microbial groups
#   48 subplot trajectories / group
#   11 annual observations (2010–2020) / trajectory
#
# MAIN OUTPUT:
#   Taylor_revision_pTTPL/pTTPL_all_groups_plot_level_full_precision.csv
#   Taylor_revision_pTTPL/pTTPL_all_groups_plot_level_full_precision.rds
#
# IMPORTANT:
# - The MPD-weighted abundance matrices are used AS PROVIDED.
# - This script does not renormalize or otherwise redefine the pTTPL input.
# - Raw numeric coefficients are never passed through format()/round() before
#   being saved or used downstream.
# - Exact R double values are preserved in the RDS output.
# - CSV numeric fields are written with up to 17 significant digits.
# ==============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------------------------
# Run this script from the repository root, or set TPL_PROJECT_ROOT to the
# project directory that contains the data folders documented in README.md.
PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

ENV_FILE <- "./data/prepared/environment_metadata.csv"

INPUT_ROOT <- "./results/03_pTTPL_MPDik"
OUTPUT_ROOT <- "./results/04_pTTPL"

GROUPS <- c(
  "Archaea", "Bacteria", "Fungi", "Protists", "Nematoda"
)

EXPECTED_YEARS <- 2010:2020
MIN_YEARS_FOR_TPL <- 5

# If FALSE, the script will stop when annual columns cannot be aligned either
# by Sample_name or by a year embedded in the column name.
# Keep FALSE for the strictest reproducibility.
ALLOW_POSITIONAL_YEAR_FALLBACK <- FALSE

# Reviewer-audit outputs. Taxon tables can be large for Bacteria.
SAVE_TAXON_TABLES <- TRUE

# A small diagnostic PDF is optional and does not affect numerical results.
SAVE_EXAMPLE_FIT_PLOTS <- TRUE
EXAMPLE_GROUP <- "Bacteria"
EXAMPLE_PLOT <- "N03"

dir.create(
  OUTPUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)

if (SAVE_TAXON_TABLES) {
  dir.create(
    file.path(OUTPUT_ROOT, "taxon_level_tables"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (SAVE_EXAMPLE_FIT_PLOTS) {
  dir.create(
    file.path(OUTPUT_ROOT, "fit_diagnostics"),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. OPTIONAL PLOTTING PACKAGE
# ------------------------------------------------------------------------------

HAS_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)

if (SAVE_EXAMPLE_FIT_PLOTS && !HAS_GGPLOT2) {
  warning(
    "ggplot2 is unavailable; example fit PDFs will be skipped."
  )
  SAVE_EXAMPLE_FIT_PLOTS <- FALSE
}

# ------------------------------------------------------------------------------
# 2. FULL-PRECISION OUTPUT HELPERS
# ------------------------------------------------------------------------------

as_full_precision_character <- function(x) {
  if (is.numeric(x)) {
    out <- sprintf("%.17g", x)
    out[is.na(x)] <- NA_character_
    out
  } else if (is.factor(x)) {
    as.character(x)
  } else {
    x
  }
}

write_csv_full_precision <- function(x, file) {

  y <- as.data.frame(
    lapply(x, as_full_precision_character),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  write.table(
    y,
    file = file,
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA",
    qmethod = "double"
  )
}

# ------------------------------------------------------------------------------
# 3. READ AND STANDARDIZE EXPERIMENTAL METADATA
# ------------------------------------------------------------------------------

env <- read.csv(
  ENV_FILE,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_env <- c(
  "Sample_name", "plot", "block", "Number", "year",
  "warm", "precip", "clip", "Treatment"
)

missing_env <- setdiff(required_env, names(env))

if (length(missing_env) > 0) {
  stop(
    "Environmental metadata are missing required columns: ",
    paste(missing_env, collapse = ", ")
  )
}

if (anyDuplicated(env$Sample_name) > 0) {
  stop("Duplicated Sample_name values were found in ENV_FILE.")
}

env$Sample_name <- as.character(env$Sample_name)
env$plot <- as.character(env$plot)
env$block <- as.character(env$block)
env$Number <- as.character(env$Number)
env$year <- as.numeric(env$year)
env$warm <- as.numeric(env$warm)
env$precip <- as.numeric(env$precip)
env$clip <- as.numeric(env$clip)

if (any(!is.finite(env$year))) {
  stop("Non-finite year values were found in ENV_FILE.")
}

# Standardize precipitation to the manuscript's continuous gradient:
#   -0.5 = half precipitation
#    0   = ambient precipitation
#   +1   = doubled precipitation
standardize_precipitation <- function(x) {

  x <- as.numeric(x)
  u <- sort(unique(x[is.finite(x)]))

  # Historical raw coding contains 0.5 and/or 2.
  if (any(abs(u - 0.5) < 1e-12) ||
      any(abs(u - 2.0) < 1e-12)) {

    out <- ifelse(
      abs(x - 0.5) < 1e-12, -0.5,
      ifelse(
        abs(x - 1.0) < 1e-12, 0,
        ifelse(
          abs(x - 2.0) < 1e-12, 1,
          NA_real_
        )
      )
    )

  } else if (
    all(u %in% c(-0.5, 0, 1))
  ) {

    out <- x

  } else {

    stop(
      "Unrecognized precipitation coding. Observed finite values: ",
      paste(u, collapse = ", ")
    )
  }

  if (any(is.na(out) & is.finite(x))) {
    stop("At least one precipitation value could not be standardized.")
  }

  out
}

env$Warming <- env$warm
env$Precipitation <- standardize_precipitation(env$precip)
env$Clipping <- env$clip

if (!all(env$Warming %in% c(0, 1))) {
  stop("Warming must be coded 0/1.")
}

if (!all(env$Clipping %in% c(0, 1))) {
  stop("Clipping must be coded 0/1.")
}

make_treatment_label <- function(w, p, c) {

  out <- rep(NA_character_, length(w))

  out[w == 0 & p ==  0   & c == 0] <- "Ctrl"
  out[w == 1 & p ==  0   & c == 0] <- "W"
  out[w == 0 & p == -0.5 & c == 0] <- "H"
  out[w == 0 & p ==  1   & c == 0] <- "D"
  out[w == 0 & p ==  0   & c == 1] <- "C"
  out[w == 1 & p == -0.5 & c == 0] <- "WH"
  out[w == 1 & p ==  1   & c == 0] <- "WD"
  out[w == 1 & p ==  0   & c == 1] <- "WC"
  out[w == 0 & p == -0.5 & c == 1] <- "HC"
  out[w == 0 & p ==  1   & c == 1] <- "DC"
  out[w == 1 & p == -0.5 & c == 1] <- "WHC"
  out[w == 1 & p ==  1   & c == 1] <- "WDC"

  out
}

env$Treatment_standard <- make_treatment_label(
  env$Warming,
  env$Precipitation,
  env$Clipping
)

if (any(is.na(env$Treatment_standard))) {
  stop("At least one Warming/Precipitation/Clipping combination is invalid.")
}

env$WholePlotID <- interaction(
  env$block,
  env$Number,
  drop = TRUE,
  sep = ":"
)

# Only these subplot IDs are accepted. Do not scan arbitrary files in folders.
PLOT_LIST <- env$plot[!duplicated(env$plot)]

message(
  "Metadata-defined subplot trajectories: n=",
  length(PLOT_LIST)
)

if (length(PLOT_LIST) != 48) {
  warning(
    "Expected 48 subplot trajectories but metadata defines ",
    length(PLOT_LIST), "."
  )
}

# ------------------------------------------------------------------------------
# 4. MATRIX I/O AND ANNUAL ALIGNMENT
# ------------------------------------------------------------------------------

read_weighted_matrix <- function(file) {

  x <- read.csv(
    file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("Empty matrix: ", file)
  }

  original_na <- sum(is.na(x))

  x[] <- lapply(
    x,
    function(z) suppressWarnings(as.numeric(z))
  )

  mat <- as.matrix(x)
  storage.mode(mat) <- "double"

  if (sum(is.na(mat)) > original_na) {
    stop(
      "Non-numeric values were introduced as NA while reading: ",
      file
    )
  }

  if (anyNA(mat)) {
    stop("NA values are present in pTTPL input: ", file)
  }

  if (any(!is.finite(mat))) {
    stop("Non-finite values are present in pTTPL input: ", file)
  }

  if (any(mat < -1e-14)) {
    stop(
      "Negative MPD-weighted abundance values were found in: ",
      file
    )
  }

  # Protect against tiny numerical negatives only.
  mat[mat < 0 & abs(mat) <= 1e-14] <- 0

  if (is.null(rownames(mat)) ||
      anyDuplicated(rownames(mat)) > 0) {
    stop(
      "Taxon row names are missing or duplicated in: ",
      file
    )
  }

  mat
}

extract_year_from_name <- function(x) {

  x <- as.character(x)

  m <- regexpr(
    "20(10|11|12|13|14|15|16|17|18|19|20)",
    x,
    perl = TRUE
  )

  out <- rep(NA_real_, length(x))
  ok <- m > 0

  out[ok] <- as.numeric(
    regmatches(x, m)[ok]
  )

  out
}

align_matrix_to_metadata <- function(
  mat,
  plot_id,
  env,
  allow_positional = FALSE
) {

  meta <- env[
    env$plot == plot_id,
    ,
    drop = FALSE
  ]

  if (nrow(meta) < MIN_YEARS_FOR_TPL) {
    stop(
      plot_id,
      ": metadata contains only ",
      nrow(meta),
      " annual records."
    )
  }

  # Check treatment / whole-plot metadata are constant across the time series.
  constant_cols <- c(
    "block", "Number", "Warming",
    "Precipitation", "Clipping",
    "Treatment_standard", "WholePlotID"
  )

  for (nm in constant_cols) {
    if (length(unique(meta[[nm]])) != 1) {
      stop(
        plot_id,
        ": metadata column '",
        nm,
        "' is not constant across years."
      )
    }
  }

  meta <- meta[
    order(meta$year),
    ,
    drop = FALSE
  ]

  expected_samples <- meta$Sample_name

  # Mode 1: exact Sample_name alignment.
  if (all(expected_samples %in% colnames(mat))) {

    used <- expected_samples

    return(list(
      mat = mat[, used, drop = FALSE],
      meta = meta,
      alignment_mode = "Sample_name_exact",
      n_input_columns = ncol(mat),
      n_used_columns = length(used)
    ))
  }

  # Mode 2: every matrix column maps to metadata by exact Sample_name, even if
  # the file contains a subset of years.
  common <- intersect(
    colnames(mat),
    expected_samples
  )

  if (length(common) >= MIN_YEARS_FOR_TPL) {

    mm <- meta[
      match(common, meta$Sample_name),
      ,
      drop = FALSE
    ]

    ord <- order(mm$year)
    mm <- mm[ord, , drop = FALSE]
    used <- mm$Sample_name

    return(list(
      mat = mat[, used, drop = FALSE],
      meta = mm,
      alignment_mode = "Sample_name_subset",
      n_input_columns = ncol(mat),
      n_used_columns = length(used)
    ))
  }

  # Mode 3: parse years directly from matrix column names.
  yrs <- extract_year_from_name(
    colnames(mat)
  )

  if (
    sum(is.finite(yrs)) >= MIN_YEARS_FOR_TPL &&
      anyDuplicated(yrs[is.finite(yrs)]) == 0
  ) {

    keep <- is.finite(yrs) &
      yrs %in% meta$year

    if (sum(keep) >= MIN_YEARS_FOR_TPL) {

      mat2 <- mat[, keep, drop = FALSE]
      yrs2 <- yrs[keep]

      ord <- order(yrs2)
      mat2 <- mat2[, ord, drop = FALSE]
      yrs2 <- yrs2[ord]

      mm <- meta[
        match(yrs2, meta$year),
        ,
        drop = FALSE
      ]

      return(list(
        mat = mat2,
        meta = mm,
        alignment_mode = "Year_parsed_from_column_name",
        n_input_columns = ncol(mat),
        n_used_columns = ncol(mat2)
      ))
    }
  }

  # Mode 4: optional positional fallback.
  if (
    isTRUE(allow_positional) &&
      ncol(mat) == nrow(meta)
  ) {

    warning(
      plot_id,
      ": annual columns were aligned positionally to sorted metadata years. ",
      "Verify that input columns are chronological."
    )

    return(list(
      mat = mat,
      meta = meta,
      alignment_mode = "POSITIONAL_FALLBACK",
      n_input_columns = ncol(mat),
      n_used_columns = ncol(mat)
    ))
  }

  stop(
    plot_id,
    ": could not align pTTPL matrix columns to annual metadata. ",
    "Expected Sample_name matching or recognizable 2010–2020 year strings. ",
    "If the file columns are known to be chronological, set ",
    "ALLOW_POSITIONAL_YEAR_FALLBACK <- TRUE after verifying the file."
  )
}

plot_metadata_one_row <- function(meta) {

  data.frame(
    block = as.character(meta$block[1]),
    Number = as.character(meta$Number[1]),
    WholePlotID = as.character(meta$WholePlotID[1]),
    warm = as.numeric(meta$warm[1]),
    precip = as.numeric(meta$precip[1]),
    clip = as.numeric(meta$clip[1]),
    Warming = as.numeric(meta$Warming[1]),
    Precipitation = as.numeric(meta$Precipitation[1]),
    Clipping = as.numeric(meta$Clipping[1]),
    Treatment = as.character(meta$Treatment_standard[1]),
    Treatment_original = as.character(meta$Treatment[1]),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. TPL FITTING HELPERS
# ------------------------------------------------------------------------------

# Vectorized linear detrending across all weighted taxon trajectories.
# The original temporal mean is retained for the detrended TPL; only temporal
# variance is replaced by the residual variance after abundance ~ year.
detrended_variance <- function(mat, years) {

  years <- as.numeric(years)

  if (length(unique(years)) < 3) {
    stop("At least three distinct years are required for detrending.")
  }

  yr <- years - mean(years)
  X <- cbind(Intercept = 1, Year = yr)

  coef_mat <- solve(
    crossprod(X),
    crossprod(X, t(mat))
  )

  fitted_mat <- t(X %*% coef_mat)
  residual_mat <- mat - fitted_mat

  # Same n-1 denominator as base R::var().
  v_resid <- rowSums(
    residual_mat^2
  ) / (ncol(mat) - 1)

  v_resid[
    v_resid < 0 &
      abs(v_resid) < 1e-15
  ] <- 0

  v_resid
}

fit_taylor <- function(
  mean_abun,
  variance_abun,
  x_ref_log10,
  taxon_names = NULL
) {

  keep <- is.finite(mean_abun) &
    is.finite(variance_abun) &
    mean_abun > 0 &
    variance_abun > 0

  if (is.null(taxon_names)) {
    taxon_names <- seq_along(mean_abun)
  }

  if (sum(keep) < 3) {

    return(list(
      Intercept = NA_real_,
      SE_Intercept = NA_real_,
      Beta = NA_real_,
      SE_Beta = NA_real_,
      Beta_CI_low = NA_real_,
      Beta_CI_high = NA_real_,
      Cov_Intercept_Beta = NA_real_,
      Cor_Intercept_Beta = NA_real_,
      R2 = NA_real_,
      Adj_R2 = NA_real_,
      Residual_SE = NA_real_,
      df_residual = NA_real_,
      F_statistic = NA_real_,
      P_Beta = NA_real_,
      LogVar_at_ref = NA_real_,
      SE_LogVar_at_ref = NA_real_,
      LogVar_at_ref_CI_low = NA_real_,
      LogVar_at_ref_CI_high = NA_real_,
      Cov_LogVarRef_Beta = NA_real_,
      Cor_LogVarRef_Beta = NA_real_,
      n_taxa = sum(keep),
      table = data.frame()
    ))
  }

  d <- data.frame(
    Taxon = as.character(taxon_names[keep]),
    Mean = mean_abun[keep],
    Variance = variance_abun[keep],
    stringsAsFactors = FALSE
  )

  d$LogMean <- log10(d$Mean)
  d$LogVar <- log10(d$Variance)

  # Original pTTPL parameterization.
  fit_raw <- lm(
    LogVar ~ LogMean,
    data = d
  )

  s_raw <- summary(fit_raw)
  vc_raw <- vcov(fit_raw)

  # Centered parameterization at the group-specific common reference.
  d$LogMean_centered <-
    d$LogMean - x_ref_log10

  fit_ref <- lm(
    LogVar ~ LogMean_centered,
    data = d
  )

  s_ref <- summary(fit_ref)
  vc_ref <- vcov(fit_ref)

  cor_raw <- vc_raw[1, 2] /
    sqrt(
      vc_raw[1, 1] *
        vc_raw[2, 2]
    )

  cor_ref <- vc_ref[1, 2] /
    sqrt(
      vc_ref[1, 1] *
        vc_ref[2, 2]
    )

  df_res <- df.residual(fit_raw)
  crit <- qt(0.975, df = df_res)

  beta <- unname(coef(fit_raw)[2])
  se_beta <- unname(
    s_raw$coefficients[2, 2]
  )

  ref_est <- unname(coef(fit_ref)[1])
  ref_se <- unname(
    s_ref$coefficients[1, 2]
  )

  fstat <- s_raw$fstatistic
  fval <- if (
    is.null(fstat)
  ) {
    NA_real_
  } else {
    unname(fstat[1])
  }

  list(
    Intercept =
      unname(coef(fit_raw)[1]),
    SE_Intercept =
      unname(s_raw$coefficients[1, 2]),
    Beta = beta,
    SE_Beta = se_beta,
    Beta_CI_low =
      beta - crit * se_beta,
    Beta_CI_high =
      beta + crit * se_beta,
    Cov_Intercept_Beta =
      unname(vc_raw[1, 2]),
    Cor_Intercept_Beta =
      unname(cor_raw),
    R2 =
      unname(s_raw$r.squared),
    Adj_R2 =
      unname(s_raw$adj.r.squared),
    Residual_SE =
      unname(s_raw$sigma),
    df_residual =
      df_res,
    F_statistic =
      fval,
    P_Beta =
      unname(s_raw$coefficients[2, 4]),
    LogVar_at_ref =
      ref_est,
    SE_LogVar_at_ref =
      ref_se,
    LogVar_at_ref_CI_low =
      ref_est - crit * ref_se,
    LogVar_at_ref_CI_high =
      ref_est + crit * ref_se,
    Cov_LogVarRef_Beta =
      unname(vc_ref[1, 2]),
    Cor_LogVarRef_Beta =
      unname(cor_ref),
    n_taxa =
      nrow(d),
    table = d
  )
}

# ------------------------------------------------------------------------------
# 6. INPUT FILE QC AND GROUP-SPECIFIC REFERENCE LEVELS
# ------------------------------------------------------------------------------

file_qc_rows <- list()
reference_rows <- list()

for (group in GROUPS) {

  message(
    "\nReference level: ",
    group
  )

  pooled_log_means <- numeric(0)

  for (plot_id in PLOT_LIST) {

    file <- file.path(
      INPUT_ROOT,
      group,
      "UPARSE",
      paste0(
        plot_id,
        "_MPD_weighted_otu_All_abun.csv"
      )
    )

    exists_flag <- file.exists(file)

    if (!exists_flag) {

      file_qc_rows[[
        length(file_qc_rows) + 1
      ]] <- data.frame(
        Group = group,
        Plot = plot_id,
        File = file,
        Exists = FALSE,
        Alignment_mode = NA_character_,
        N_taxa = NA_integer_,
        N_input_columns = NA_integer_,
        N_used_columns = NA_integer_,
        First_year = NA_real_,
        Last_year = NA_real_,
        stringsAsFactors = FALSE
      )

      next
    }

    mat <- read_weighted_matrix(file)

    ali <- align_matrix_to_metadata(
      mat,
      plot_id,
      env,
      allow_positional =
        ALLOW_POSITIONAL_YEAR_FALLBACK
    )

    mat <- ali$mat
    meta <- ali$meta

    mu <- rowMeans(mat)
    ok <- is.finite(mu) & mu > 0

    pooled_log_means <- c(
      pooled_log_means,
      log10(mu[ok])
    )

    file_qc_rows[[
      length(file_qc_rows) + 1
    ]] <- data.frame(
      Group = group,
      Plot = plot_id,
      File = file,
      Exists = TRUE,
      Alignment_mode =
        ali$alignment_mode,
      N_taxa = nrow(mat),
      N_input_columns =
        ali$n_input_columns,
      N_used_columns =
        ali$n_used_columns,
      First_year = min(meta$year),
      Last_year = max(meta$year),
      stringsAsFactors = FALSE
    )
  }

  n_missing <- sum(
    vapply(
      PLOT_LIST,
      function(plot_id) {
        !file.exists(
          file.path(
            INPUT_ROOT,
            group,
            "UPARSE",
            paste0(
              plot_id,
              "_MPD_weighted_otu_All_abun.csv"
            )
          )
        )
      },
      logical(1)
    )
  )

  if (n_missing > 0) {

    # Save the QC accumulated so far before stopping.
    current_qc <- do.call(
      rbind,
      file_qc_rows
    )

    write_csv_full_precision(
      current_qc,
      file.path(
        OUTPUT_ROOT,
        "00_pTTPL_input_file_QC.csv"
      )
    )

    stop(
      group,
      ": missing ",
      n_missing,
      " expected pTTPL input file(s). ",
      "See Taylor_revision_pTTPL/00_pTTPL_input_file_QC.csv."
    )
  }

  if (length(pooled_log_means) == 0) {
    stop(
      "No positive finite temporal means were available for ",
      group
    )
  }

  xref <- median(
    pooled_log_means,
    na.rm = TRUE
  )

  reference_rows[[
    length(reference_rows) + 1
  ]] <- data.frame(
    Group = group,
    Log10WeightedMeanReference = xref,
    WeightedMeanReference = 10^xref,
    N_pooled_taxon_plot_means =
      length(pooled_log_means),
    stringsAsFactors = FALSE
  )
}

file_qc <- do.call(
  rbind,
  file_qc_rows
)

xref_table <- do.call(
  rbind,
  reference_rows
)

write_csv_full_precision(
  file_qc,
  file.path(
    OUTPUT_ROOT,
    "00_pTTPL_input_file_QC.csv"
  )
)

write_csv_full_precision(
  xref_table,
  file.path(
    OUTPUT_ROOT,
    "01_pTTPL_group_reference_level_full_precision.csv"
  )
)

saveRDS(
  xref_table,
  file.path(
    OUTPUT_ROOT,
    "01_pTTPL_group_reference_level_full_precision.rds"
  ),
  version = 3
)

# ------------------------------------------------------------------------------
# 7. RE-FIT ALL 48 x 5 pTTPL TRAJECTORIES
# ------------------------------------------------------------------------------

all_rows <- list()
counter <- 1

for (group in GROUPS) {

  message(
    "\n======================================"
  )
  message(
    "Fitting pTTPL: ",
    group
  )
  message(
    "======================================"
  )

  x_ref <- xref_table$
    Log10WeightedMeanReference[
      xref_table$Group == group
    ]

  group_taxon_dir <- file.path(
    OUTPUT_ROOT,
    "taxon_level_tables",
    group
  )

  if (SAVE_TAXON_TABLES) {
    dir.create(
      group_taxon_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  for (plot_id in PLOT_LIST) {

    message(
      "  ",
      group,
      " / ",
      plot_id
    )

    file <- file.path(
      INPUT_ROOT,
      group,
      "UPARSE",
      paste0(
        plot_id,
        "_MPD_weighted_otu_All_abun.csv"
      )
    )

    mat <- read_weighted_matrix(file)

    ali <- align_matrix_to_metadata(
      mat,
      plot_id,
      env,
      allow_positional =
        ALLOW_POSITIONAL_YEAR_FALLBACK
    )

    mat <- ali$mat
    meta <- ali$meta
    years <- meta$year

    mean_weighted <- rowMeans(mat)
    var_raw <- apply(
      mat,
      1,
      var
    )

    var_detrended <- detrended_variance(
      mat,
      years
    )

    tpl_raw <- fit_taylor(
      mean_weighted,
      var_raw,
      x_ref_log10 = x_ref,
      taxon_names = rownames(mat)
    )

    tpl_det <- fit_taylor(
      mean_weighted,
      var_detrended,
      x_ref_log10 = x_ref,
      taxon_names = rownames(mat)
    )

    md <- plot_metadata_one_row(meta)

    result <- cbind(
      data.frame(
        TPL_type = "pTTPL",
        Input_metric =
          "MPD_weighted_abundance",
        Group = group,
        Plot = plot_id,
        n_years = ncol(mat),
        n_taxa_total = nrow(mat),
        Alignment_mode =
          ali$alignment_mode,
        Log10WeightedMeanReference =
          x_ref,
        WeightedMeanReference =
          10^x_ref,
        stringsAsFactors = FALSE
      ),
      md,
      data.frame(
        # ------------------------------------------------------
        # RAW pTTPL
        # ------------------------------------------------------
        Intercept_raw =
          tpl_raw$Intercept,
        Alpha_raw_scale =
          10^tpl_raw$Intercept,
        SE_Intercept_raw =
          tpl_raw$SE_Intercept,
        Beta_raw =
          tpl_raw$Beta,
        SE_Beta_raw =
          tpl_raw$SE_Beta,
        Beta_CI_low_raw =
          tpl_raw$Beta_CI_low,
        Beta_CI_high_raw =
          tpl_raw$Beta_CI_high,
        Cov_Intercept_Beta_raw =
          tpl_raw$Cov_Intercept_Beta,
        Cor_Intercept_Beta_raw =
          tpl_raw$Cor_Intercept_Beta,
        R2_raw =
          tpl_raw$R2,
        Adj_R2_raw =
          tpl_raw$Adj_R2,
        Residual_SE_raw =
          tpl_raw$Residual_SE,
        df_residual_raw =
          tpl_raw$df_residual,
        F_statistic_raw =
          tpl_raw$F_statistic,
        P_Beta_raw =
          tpl_raw$P_Beta,
        n_taxa_TPL_raw =
          tpl_raw$n_taxa,

        # Centered raw pTTPL
        LogVar_at_ref_raw =
          tpl_raw$LogVar_at_ref,
        Var_at_ref_raw =
          10^tpl_raw$LogVar_at_ref,
        SE_LogVar_at_ref_raw =
          tpl_raw$SE_LogVar_at_ref,
        LogVar_at_ref_CI_low_raw =
          tpl_raw$LogVar_at_ref_CI_low,
        LogVar_at_ref_CI_high_raw =
          tpl_raw$LogVar_at_ref_CI_high,
        Cov_LogVarRef_Beta_raw =
          tpl_raw$Cov_LogVarRef_Beta,
        Cor_LogVarRef_Beta_raw =
          tpl_raw$Cor_LogVarRef_Beta,

        # ------------------------------------------------------
        # DETRENDED pTTPL sensitivity
        # ------------------------------------------------------
        Intercept_detrended =
          tpl_det$Intercept,
        Alpha_detrended_scale =
          10^tpl_det$Intercept,
        SE_Intercept_detrended =
          tpl_det$SE_Intercept,
        Beta_detrended =
          tpl_det$Beta,
        SE_Beta_detrended =
          tpl_det$SE_Beta,
        Beta_CI_low_detrended =
          tpl_det$Beta_CI_low,
        Beta_CI_high_detrended =
          tpl_det$Beta_CI_high,
        Cov_Intercept_Beta_detrended =
          tpl_det$Cov_Intercept_Beta,
        Cor_Intercept_Beta_detrended =
          tpl_det$Cor_Intercept_Beta,
        R2_detrended =
          tpl_det$R2,
        Adj_R2_detrended =
          tpl_det$Adj_R2,
        Residual_SE_detrended =
          tpl_det$Residual_SE,
        df_residual_detrended =
          tpl_det$df_residual,
        F_statistic_detrended =
          tpl_det$F_statistic,
        P_Beta_detrended =
          tpl_det$P_Beta,
        n_taxa_TPL_detrended =
          tpl_det$n_taxa,

        # Centered detrended pTTPL
        LogVar_at_ref_detrended =
          tpl_det$LogVar_at_ref,
        Var_at_ref_detrended =
          10^tpl_det$LogVar_at_ref,
        SE_LogVar_at_ref_detrended =
          tpl_det$SE_LogVar_at_ref,
        LogVar_at_ref_CI_low_detrended =
          tpl_det$LogVar_at_ref_CI_low,
        LogVar_at_ref_CI_high_detrended =
          tpl_det$LogVar_at_ref_CI_high,
        Cov_LogVarRef_Beta_detrended =
          tpl_det$Cov_LogVarRef_Beta,
        Cor_LogVarRef_Beta_detrended =
          tpl_det$Cor_LogVarRef_Beta,

        stringsAsFactors = FALSE
      )
    )

    all_rows[[counter]] <- result
    counter <- counter + 1

    # ----------------------------------------------------------
    # Optional per-plot taxon-level audit table.
    # ----------------------------------------------------------
    if (SAVE_TAXON_TABLES) {

      taxon_audit <- data.frame(
        Taxon = rownames(mat),
        Temporal_mean_weighted_abundance =
          mean_weighted,
        Temporal_variance_raw =
          var_raw,
        Temporal_variance_detrended =
          var_detrended,
        Log10_mean_weighted_abundance =
          ifelse(
            mean_weighted > 0,
            log10(mean_weighted),
            NA_real_
          ),
        Log10_variance_raw =
          ifelse(
            var_raw > 0,
            log10(var_raw),
            NA_real_
          ),
        Log10_variance_detrended =
          ifelse(
            var_detrended > 0,
            log10(var_detrended),
            NA_real_
          ),
        Used_in_raw_pTTPL =
          is.finite(mean_weighted) &
          mean_weighted > 0 &
          is.finite(var_raw) &
          var_raw > 0,
        Used_in_detrended_pTTPL =
          is.finite(mean_weighted) &
          mean_weighted > 0 &
          is.finite(var_detrended) &
          var_detrended > 0,
        stringsAsFactors = FALSE
      )

      write_csv_full_precision(
        taxon_audit,
        file.path(
          group_taxon_dir,
          paste0(
            plot_id,
            "_pTTPL_taxon_mean_variance_full_precision.csv"
          )
        )
      )
    }

    # ----------------------------------------------------------
    # Optional example fit plot.
    # ----------------------------------------------------------
    if (
      SAVE_EXAMPLE_FIT_PLOTS &&
      identical(group, EXAMPLE_GROUP) &&
      identical(plot_id, EXAMPLE_PLOT) &&
      nrow(tpl_raw$table) > 0
    ) {

      p <- ggplot2::ggplot(
        tpl_raw$table,
        ggplot2::aes(
          x = LogMean,
          y = LogVar
        )
      ) +
        ggplot2::geom_point(
          size = 1,
          alpha = 0.55
        ) +
        ggplot2::geom_smooth(
          method = "lm",
          formula = y ~ x,
          se = TRUE,
          linewidth = 0.8
        ) +
        ggplot2::theme_bw() +
        ggplot2::labs(
          title = paste(
            EXAMPLE_GROUP,
            EXAMPLE_PLOT,
            "pTTPL"
          ),
          x =
            "log10 temporal mean MPD-weighted abundance",
          y =
            "log10 temporal variance"
        )

      ggplot2::ggsave(
        file.path(
          OUTPUT_ROOT,
          "fit_diagnostics",
          paste0(
            EXAMPLE_GROUP,
            "_",
            EXAMPLE_PLOT,
            "_pTTPL_example.pdf"
          )
        ),
        p,
        width = 5,
        height = 4.5
      )
    }
  }
}

all_metrics <- do.call(
  rbind,
  all_rows
)

# Fixed manuscript-oriented group ordering.
all_metrics$Group <- factor(
  all_metrics$Group,
  levels = GROUPS
)

all_metrics$Treatment <- factor(
  all_metrics$Treatment,
  levels = c(
    "Ctrl", "W", "H", "D", "C",
    "WH", "WD", "WC", "HC", "DC",
    "WHC", "WDC"
  )
)

all_metrics <- all_metrics[
  order(
    all_metrics$Group,
    match(
      all_metrics$Plot,
      PLOT_LIST
    )
  ),
  ,
  drop = FALSE
]

all_metrics$Group <- as.character(
  all_metrics$Group
)
all_metrics$Treatment <- as.character(
  all_metrics$Treatment
)

# ------------------------------------------------------------------------------
# 8. SAVE EXACT / FULL-PRECISION RESULTS
# ------------------------------------------------------------------------------

saveRDS(
  all_metrics,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_all_groups_plot_level_full_precision.rds"
  ),
  version = 3
)

write_csv_full_precision(
  all_metrics,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_all_groups_plot_level_full_precision.csv"
  )
)

# One file per microbial group.
for (group in GROUPS) {

  dg <- all_metrics[
    all_metrics$Group == group,
    ,
    drop = FALSE
  ]

  write_csv_full_precision(
    dg,
    file.path(
      OUTPUT_ROOT,
      paste0(
        "pTTPL_",
        group,
        "_plot_level_full_precision.csv"
      )
    )
  )
}

# ------------------------------------------------------------------------------
# 9. REPLACEMENT SUPPLEMENTARY TABLE 2 INPUT
# ------------------------------------------------------------------------------

supp_long <- all_metrics[
  ,
  c(
    "Treatment", "Plot", "Group",
    "Intercept_raw",
    "Alpha_raw_scale",
    "Beta_raw",
    "SE_Beta_raw",
    "Beta_CI_low_raw",
    "Beta_CI_high_raw",
    "R2_raw",
    "Adj_R2_raw",
    "P_Beta_raw",
    "n_taxa_TPL_raw"
  ),
  drop = FALSE
]

names(supp_long)[
  names(supp_long) == "Beta_raw"
] <- "Beta_p"

names(supp_long)[
  names(supp_long) == "SE_Beta_raw"
] <- "SE_Beta_p"

write_csv_full_precision(
  supp_long,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_SuppTable2_long_full_precision.csv"
  )
)

# Wide table: one row per Treatment x Plot, with a block of exact pTTPL
# parameters for each microbial group.
wide_parts <- lapply(
  GROUPS,
  function(group) {

    d <- supp_long[
      supp_long$Group == group,
      ,
      drop = FALSE
    ]

    keep <- c(
      "Treatment", "Plot",
      "Intercept_raw",
      "Alpha_raw_scale",
      "Beta_p",
      "SE_Beta_p",
      "Beta_CI_low_raw",
      "Beta_CI_high_raw",
      "R2_raw",
      "P_Beta_raw"
    )

    d <- d[, keep, drop = FALSE]

    parameter_cols <- setdiff(
      names(d),
      c("Treatment", "Plot")
    )

    names(d)[
      match(parameter_cols, names(d))
    ] <- paste0(
      group,
      "_",
      parameter_cols
    )

    d
  }
)

supp_wide <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = c("Treatment", "Plot"),
      all = TRUE,
      sort = FALSE
    )
  },
  wide_parts
)

supp_wide$Treatment_order <- match(
  supp_wide$Treatment,
  c(
    "Ctrl", "W", "H", "D", "C",
    "WH", "WD", "WC", "HC", "DC",
    "WHC", "WDC"
  )
)

supp_wide$Plot_order <- match(
  supp_wide$Plot,
  PLOT_LIST
)

supp_wide <- supp_wide[
  order(
    supp_wide$Plot_order
  ),
  ,
  drop = FALSE
]

supp_wide$Treatment_order <- NULL
supp_wide$Plot_order <- NULL

write_csv_full_precision(
  supp_wide,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_SuppTable2_wide_full_precision.csv"
  )
)

# A manuscript-friendly rounded VIEW only.
# Do not use this rounded file as downstream statistical input.
supp_view <- supp_wide

num_cols <- vapply(
  supp_view,
  is.numeric,
  logical(1)
)

supp_view[num_cols] <- lapply(
  supp_view[num_cols],
  function(x) round(x, 6)
)

write.csv(
  supp_view,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_SuppTable2_manuscript_view_6dp.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------------------
# 10. FIT QC SUMMARY
# ------------------------------------------------------------------------------

fit_qc <- aggregate(
  cbind(
    R2_raw,
    R2_detrended,
    SE_Beta_raw,
    SE_Beta_detrended,
    n_taxa_TPL_raw,
    n_taxa_TPL_detrended
  ) ~ Group,
  data = all_metrics,
  FUN = function(x) {
    c(
      min = min(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
)

saveRDS(
  fit_qc,
  file.path(
    OUTPUT_ROOT,
    "pTTPL_fit_QC_summary.rds"
  ),
  version = 3
)

# Simpler long-form QC table.
qc_rows <- list()
qq <- 1

for (group in GROUPS) {

  d <- all_metrics[
    all_metrics$Group == group,
    ,
    drop = FALSE
  ]

  qc_rows[[qq]] <- data.frame(
    Group = group,
    N_plots = nrow(d),
    N_blocks = length(
      unique(d$block)
    ),
    N_wholeplots = length(
      unique(d$WholePlotID)
    ),
    R2_raw_min =
      min(d$R2_raw, na.rm = TRUE),
    R2_raw_median =
      median(d$R2_raw, na.rm = TRUE),
    R2_raw_max =
      max(d$R2_raw, na.rm = TRUE),
    Beta_raw_min =
      min(d$Beta_raw, na.rm = TRUE),
    Beta_raw_median =
      median(d$Beta_raw, na.rm = TRUE),
    Beta_raw_max =
      max(d$Beta_raw, na.rm = TRUE),
    SE_Beta_raw_min =
      min(d$SE_Beta_raw, na.rm = TRUE),
    SE_Beta_raw_median =
      median(d$SE_Beta_raw, na.rm = TRUE),
    SE_Beta_raw_max =
      max(d$SE_Beta_raw, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  qq <- qq + 1
}

qc_simple <- do.call(
  rbind,
  qc_rows
)

write_csv_full_precision(
  qc_simple,
  file.path(
    OUTPUT_ROOT,
    "04_pTTPL_fit_QC_summary.csv"
  )
)

# ------------------------------------------------------------------------------
# 11. SESSION INFO / README
# ------------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_ROOT,
    "sessionInfo_pTTPL_fit.txt"
  )
)

readme <- c(
  "NWS pTTPL full-precision refit",
  "",
  "PRIMARY NUMERICAL INPUT FOR DOWNSTREAM STATISTICS:",
  "pTTPL_all_groups_plot_level_full_precision.rds",
  "pTTPL_all_groups_plot_level_full_precision.csv",
  "",
  "The RDS preserves exact R double values.",
  "The CSV writes numeric columns with up to 17 significant digits.",
  "",
  "Do NOT use pTTPL_SuppTable2_manuscript_view_6dp.csv for statistics.",
  "That file is only a human-readable rounded view.",
  "",
  "Core fitted quantities:",
  "Beta_raw, SE_Beta_raw, beta 95% CI, intercept and SE, R2, adjusted R2,",
  "slope P, centered reference-level variance, and full covariance terms.",
  "",
  "Detrended sensitivity:",
  "Weighted abundance was linearly detrended within each taxon-by-subplot",
  "trajectory; residual temporal variance was used while retaining the original",
  "11-year temporal mean for the TPL x axis.",
  "",
  "NEXT STEP:",
  "Run NWS_pTTPL_02_complete_revision_statistics.R."
)

writeLines(
  readme,
  file.path(
    OUTPUT_ROOT,
    "README_pTTPL_fit.txt"
  )
)

message(
  "\n============================================================"
)
message(
  "pTTPL full-precision refit completed."
)
message(
  "Primary exact output:\n  ",
  file.path(
    OUTPUT_ROOT,
    "pTTPL_all_groups_plot_level_full_precision.rds"
  )
)
message(
  "Primary CSV output:\n  ",
  file.path(
    OUTPUT_ROOT,
    "pTTPL_all_groups_plot_level_full_precision.csv"
  )
)
message(
  "============================================================\n"
)
