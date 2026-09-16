# ==============================================================================
# NWS pTTPL SUBGROUP WORKFLOW — STEP 1
# Full-precision phylum / trophic-level pTTPL fitting, 2010–2020
#
# Purpose
#   1) Refit subplot-level pTTPL beta within phyla / trophic groups using the
#      existing MPD-weighted OTU matrices from the current pTTPL pipeline.
#   2) Preserve the MPD-weighted values AS PROVIDED: no renormalization.
#   3) Use taxonomy / trophic files only to define subgroup OTU membership.
#   4) Apply the minimum-taxon requirement to taxa that ACTUALLY enter the TPL
#      regression (finite weighted mean > 0 and temporal variance > 0).
#   5) Keep TPL-fit reliability separate from the later Fig. 2d display filter.
#
# Weighted pTTPL input:
#   ./TPL_Phylogenetic_beta/predata/{Group}/UPARSE/
#       {Plot}_MPD_weighted_otu_All_abun.csv
#
# Membership input:
#   Phylum:
#   ./predata/2010to2020/{Group}/
#       {Plot}_UPARSE_otu_All_taxonomy_abun_freq.csv
#
#   Trophic:
#   ./predata/2010to2020/Protists/Trophic/
#       Trophic_{level}_{Plot}_UPARSE_otu_All.csv
#   ./predata/2010to2020/Nematoda/Trophic/
#       Trophic_{level}_{Plot}_UPARSE_otu_All.csv
#
# Main outputs
#   Taylor_revision_pTTPL_subgroup/
#     pTTPL_phylum_trophic_plot_level_full_precision.csv
#     pTTPL_phylum_trophic_plot_level_full_precision.rds
#     pTTPL_phylum_trophic_input_QC.csv
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

WEIGHTED_INPUT_ROOT <- "./results/03_pTTPL_MPDik"
MEMBERSHIP_ROOT <- "./data/prepared/taxonomy"
OUTPUT_ROOT <- "./results/08_pTTPL_subgroup"

GROUPS <- c("Archaea", "Bacteria", "Fungi", "Protists", "Nematoda")
PHYLUM_GROUPS <- c("Archaea", "Bacteria", "Fungi", "Protists")

TROPHIC_LEVELS <- list(
  Nematoda = c("b", "f", "h", "o", "p")
)

MIN_YEARS <- 5L

# Reliability thresholds for a SINGLE subplot-level TPL fit.
MIN_TAXA_PHYLUM  <- 15L
MIN_TAXA_TROPHIC <- 6L

# Strict annual alignment follows the current pTTPL workflow.
ALLOW_POSITIONAL_YEAR_FALLBACK <- FALSE

SAVE_SUBGROUP_TAXON_TABLES <- FALSE

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)
if (SAVE_SUBGROUP_TAXON_TABLES) {
  dir.create(file.path(OUTPUT_ROOT, "taxon_tables"),
             recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# 1. FULL-PRECISION OUTPUT HELPERS
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
    y, file = file, sep = ",",
    row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = "NA", qmethod = "double"
  )
}

# ------------------------------------------------------------------------------
# 2. METADATA AND EXPERIMENTAL DESIGN
# ------------------------------------------------------------------------------

env <- read.csv(
  ENV_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)

required_env <- c(
  "Sample_name", "block", "plot", "Number", "year",
  "warm", "precip", "clip", "Treatment"
)

miss <- setdiff(required_env, names(env))
if (length(miss) > 0) {
  stop("ENV_FILE is missing required columns: ", paste(miss, collapse = ", "))
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

env$WarmF <- factor(env$warm, levels = c(0, 1),
                    labels = c("Unwarmed", "Warmed"))

env$PrecipF <- factor(
  env$precip,
  levels = c(1, 0.5, 2),
  labels = c("Ambient", "Half", "Double")
)

env$ClipF <- factor(env$clip, levels = c(0, 1),
                    labels = c("Unclipped", "Clipped"))

env$WholePlotID <- interaction(
  env$block, env$Number, drop = TRUE, sep = ":"
)

make_treatment_label <- function(warm, precip, clip) {
  out <- rep(NA_character_, length(warm))
  out[warm == 0 & precip == 1   & clip == 0] <- "Ctrl"
  out[warm == 1 & precip == 1   & clip == 0] <- "W"
  out[warm == 0 & precip == 0.5 & clip == 0] <- "H"
  out[warm == 0 & precip == 2   & clip == 0] <- "D"
  out[warm == 0 & precip == 1   & clip == 1] <- "C"
  out[warm == 1 & precip == 0.5 & clip == 0] <- "WH"
  out[warm == 1 & precip == 2   & clip == 0] <- "WD"
  out[warm == 1 & precip == 1   & clip == 1] <- "WC"
  out[warm == 0 & precip == 0.5 & clip == 1] <- "HC"
  out[warm == 0 & precip == 2   & clip == 1] <- "DC"
  out[warm == 1 & precip == 0.5 & clip == 1] <- "WHC"
  out[warm == 1 & precip == 2   & clip == 1] <- "WDC"
  out
}

env$Treatment_standard <- make_treatment_label(
  env$warm, env$precip, env$clip
)

if (any(is.na(env$Treatment_standard))) {
  stop("Invalid warm/precip/clip combination detected in metadata.")
}

PLOT_LIST <- env$plot[!duplicated(env$plot)]
message("Metadata-defined subplot trajectories: n = ", length(PLOT_LIST))

# ------------------------------------------------------------------------------
# 3. MATRIX I/O AND ANNUAL ALIGNMENT
# ------------------------------------------------------------------------------

read_weighted_matrix <- function(file) {
  x <- read.csv(
    file, header = TRUE, row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  if (nrow(x) == 0 || ncol(x) == 0) stop("Empty matrix: ", file)

  original_na <- sum(is.na(x))
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  mat <- as.matrix(x)
  storage.mode(mat) <- "double"

  if (sum(is.na(mat)) > original_na) {
    stop("Non-numeric values were introduced as NA while reading: ", file)
  }
  if (anyNA(mat)) stop("NA values are present in pTTPL input: ", file)
  if (any(!is.finite(mat))) stop("Non-finite values are present in pTTPL input: ", file)
  if (any(mat < -1e-14)) stop("Negative weighted values found in: ", file)

  mat[mat < 0 & abs(mat) <= 1e-14] <- 0

  if (is.null(rownames(mat)) || anyDuplicated(rownames(mat)) > 0) {
    stop("Taxon row names are missing or duplicated in: ", file)
  }

  mat
}

extract_year_from_name <- function(x) {
  x <- as.character(x)
  m <- regexpr("20(10|11|12|13|14|15|16|17|18|19|20)", x, perl = TRUE)
  out <- rep(NA_real_, length(x))
  ok <- m > 0
  out[ok] <- as.numeric(regmatches(x, m)[ok])
  out
}

align_matrix_to_metadata <- function(
  mat, plot_id, env, allow_positional = FALSE
) {
  meta <- env[env$plot == plot_id, , drop = FALSE]

  if (nrow(meta) < MIN_YEARS) {
    stop(plot_id, ": fewer than ", MIN_YEARS, " metadata years.")
  }

  const_cols <- c(
    "block", "Number", "warm", "precip", "clip",
    "Treatment_standard", "WarmF", "PrecipF", "ClipF", "WholePlotID"
  )

  for (nm in const_cols) {
    if (length(unique(meta[[nm]])) != 1) {
      stop(plot_id, ": metadata column ", nm, " varies among years.")
    }
  }

  meta <- meta[order(meta$year), , drop = FALSE]
  expected <- meta$Sample_name

  # Mode 1: exact Sample_name alignment
  if (all(expected %in% colnames(mat))) {
    return(list(
      mat = mat[, expected, drop = FALSE],
      meta = meta,
      alignment_mode = "Sample_name_exact"
    ))
  }

  # Mode 2: Sample_name subset
  common <- intersect(colnames(mat), expected)
  if (length(common) >= MIN_YEARS) {
    mm <- meta[match(common, meta$Sample_name), , drop = FALSE]
    ord <- order(mm$year)
    mm <- mm[ord, , drop = FALSE]

    return(list(
      mat = mat[, mm$Sample_name, drop = FALSE],
      meta = mm,
      alignment_mode = "Sample_name_subset"
    ))
  }

  # Mode 3: parse years from matrix column names
  yrs <- extract_year_from_name(colnames(mat))
  keep <- is.finite(yrs) & yrs %in% meta$year

  if (sum(keep) >= MIN_YEARS &&
      anyDuplicated(yrs[keep]) == 0) {

    mat2 <- mat[, keep, drop = FALSE]
    yrs2 <- yrs[keep]
    ord <- order(yrs2)
    mat2 <- mat2[, ord, drop = FALSE]
    yrs2 <- yrs2[ord]

    mm <- meta[match(yrs2, meta$year), , drop = FALSE]

    return(list(
      mat = mat2,
      meta = mm,
      alignment_mode = "Year_parsed_from_column_name"
    ))
  }

  # Mode 4: optional positional fallback
  if (isTRUE(allow_positional) && ncol(mat) == nrow(meta)) {
    warning(plot_id, ": using positional year alignment after verification.")
    return(list(
      mat = mat,
      meta = meta,
      alignment_mode = "POSITIONAL_FALLBACK"
    ))
  }

  stop(
    plot_id,
    ": could not align pTTPL annual columns to metadata. ",
    "Keep ALLOW_POSITIONAL_YEAR_FALLBACK=FALSE unless chronology is verified."
  )
}

plot_metadata_one_row <- function(meta) {
  data.frame(
    block = as.character(meta$block[1]),
    Plot = as.character(meta$plot[1]),
    Number = as.character(meta$Number[1]),
    WholePlotID = as.character(meta$WholePlotID[1]),
    warm = as.numeric(meta$warm[1]),
    precip = as.numeric(meta$precip[1]),
    clip = as.numeric(meta$clip[1]),
    WarmF = as.character(meta$WarmF[1]),
    PrecipF = as.character(meta$PrecipF[1]),
    ClipF = as.character(meta$ClipF[1]),
    Treatment = as.character(meta$Treatment_standard[1]),
    stringsAsFactors = FALSE
  )
}

# Global subgroup-membership tables are stored once per soil-biotic group in
# ./data/prepared/taxonomy. Membership does not vary among plots.
TAXONOMY_FILES <- c(
  Archaea  = file.path("./data/prepared/taxonomy", "Archaea_taxonomy.csv"),
  Bacteria = file.path("./data/prepared/taxonomy", "Bacteria_taxonomy.csv"),
  Fungi    = file.path("./data/prepared/taxonomy", "Fungi_taxonomy.csv"),
  Protists = file.path("./data/prepared/taxonomy", "Protists_taxonomy.csv"),
  Nematoda = file.path("./data/prepared/taxonomy", "Nematoda_taxonomy.csv")
)

read_taxonomy <- function(group, plot_id = NULL) {
  file <- unname(TAXONOMY_FILES[group])
  if (length(file) != 1L || is.na(file) || !file.exists(file)) return(NULL)

  tx <- read.csv(
    file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!"TaxonID" %in% names(tx)) stop("Taxonomy file lacks TaxonID: ", file)
  ids <- as.character(tx$TaxonID)
  tx$TaxonID <- NULL
  rownames(tx) <- ids
  if (anyDuplicated(rownames(tx)) > 0) {
    stop("Duplicated TaxonID values in: ", file)
  }
  tx
}

read_trophic_membership <- function(group, trophic) {
  # The uploaded primary data contain trophic assignments for Nematoda.
  # Protist trophic assignments were not included in the released input set, so
  # the canonical subgroup workflow uses Protist phyla and Nematoda trophic groups.
  if (group != "Nematoda") return(character(0))
  tx <- read_taxonomy("Nematoda")
  if (is.null(tx) || !"trophic_groups" %in% names(tx)) return(character(0))
  rownames(tx)[as.character(tx$trophic_groups) == trophic]
}

# ------------------------------------------------------------------------------
# 4. TPL FITTING HELPERS
# ------------------------------------------------------------------------------

detrended_variance <- function(mat, years) {
  years <- as.numeric(years)
  if (length(unique(years)) < 3) {
    return(rep(NA_real_, nrow(mat)))
  }

  yr <- years - mean(years)
  X <- cbind(Intercept = 1, Year = yr)
  coef_mat <- solve(crossprod(X), crossprod(X, t(mat)))
  fitted_mat <- t(X %*% coef_mat)
  residual_mat <- mat - fitted_mat

  v_resid <- rowSums(residual_mat^2) / (ncol(mat) - 1)
  v_resid[v_resid < 0 & abs(v_resid) < 1e-15] <- 0
  v_resid
}

fit_one_taylor <- function(mean_abun, variance_abun, min_taxa) {
  keep <- is.finite(mean_abun) &
    is.finite(variance_abun) &
    mean_abun > 0 &
    variance_abun > 0

  n_used <- sum(keep)

  empty <- list(
    Fit_ok = FALSE,
    Intercept = NA_real_,
    SE_Intercept = NA_real_,
    Beta = NA_real_,
    SE_Beta = NA_real_,
    Beta_CI_low = NA_real_,
    Beta_CI_high = NA_real_,
    R2 = NA_real_,
    Adj_R2 = NA_real_,
    P_Beta = NA_real_,
    n_taxa_used = n_used,
    table = data.frame()
  )

  if (n_used < min_taxa) return(empty)

  d <- data.frame(
    Mean = mean_abun[keep],
    Variance = variance_abun[keep]
  )
  d$LogMean <- log10(d$Mean)
  d$LogVar <- log10(d$Variance)

  fit <- lm(LogVar ~ LogMean, data = d)
  sm <- summary(fit)
  ci <- suppressMessages(confint(fit, "LogMean", level = 0.95))

  list(
    Fit_ok = TRUE,
    Intercept = unname(coef(fit)[1]),
    SE_Intercept = unname(sm$coefficients[1, 2]),
    Beta = unname(coef(fit)[2]),
    SE_Beta = unname(sm$coefficients[2, 2]),
    Beta_CI_low = unname(ci[1]),
    Beta_CI_high = unname(ci[2]),
    R2 = unname(sm$r.squared),
    Adj_R2 = unname(sm$adj.r.squared),
    P_Beta = unname(sm$coefficients[2, 4]),
    n_taxa_used = n_used,
    table = d
  )
}

fit_subgroup <- function(mat_sub, years, min_taxa) {
  mu <- rowMeans(mat_sub)
  vv_raw <- apply(mat_sub, 1, var)
  vv_det <- detrended_variance(mat_sub, years)

  raw <- fit_one_taylor(mu, vv_raw, min_taxa)
  det <- fit_one_taylor(mu, vv_det, min_taxa)

  list(raw = raw, detrended = det)
}

make_result_row <- function(
  group, level_type, subgroup, meta_one,
  n_membership, n_matched, fit_obj, min_taxa,
  alignment_mode
) {
  data.frame(
    TPL_type = "pTTPL",
    Group = group,
    Level_type = level_type,
    Subgroup = subgroup,
    Plot = meta_one$Plot,
    block = meta_one$block,
    Number = meta_one$Number,
    WholePlotID = meta_one$WholePlotID,
    warm = meta_one$warm,
    precip = meta_one$precip,
    clip = meta_one$clip,
    WarmF = meta_one$WarmF,
    PrecipF = meta_one$PrecipF,
    ClipF = meta_one$ClipF,
    Treatment = meta_one$Treatment,
    Alignment_mode = alignment_mode,
    Min_taxa_required = min_taxa,
    n_taxa_membership = n_membership,
    n_taxa_matched = n_matched,

    Fit_ok_raw = fit_obj$raw$Fit_ok,
    n_taxa_used_raw = fit_obj$raw$n_taxa_used,
    Intercept_raw = fit_obj$raw$Intercept,
    SE_Intercept_raw = fit_obj$raw$SE_Intercept,
    Beta_raw = fit_obj$raw$Beta,
    SE_Beta_raw = fit_obj$raw$SE_Beta,
    Beta_raw_CI_low = fit_obj$raw$Beta_CI_low,
    Beta_raw_CI_high = fit_obj$raw$Beta_CI_high,
    R2_raw = fit_obj$raw$R2,
    Adj_R2_raw = fit_obj$raw$Adj_R2,
    P_Beta_raw = fit_obj$raw$P_Beta,

    Fit_ok_detrended = fit_obj$detrended$Fit_ok,
    n_taxa_used_detrended = fit_obj$detrended$n_taxa_used,
    Intercept_detrended = fit_obj$detrended$Intercept,
    SE_Intercept_detrended = fit_obj$detrended$SE_Intercept,
    Beta_detrended = fit_obj$detrended$Beta,
    SE_Beta_detrended = fit_obj$detrended$SE_Beta,
    Beta_detrended_CI_low = fit_obj$detrended$Beta_CI_low,
    Beta_detrended_CI_high = fit_obj$detrended$Beta_CI_high,
    R2_detrended = fit_obj$detrended$R2,
    Adj_R2_detrended = fit_obj$detrended$Adj_R2,
    P_Beta_detrended = fit_obj$detrended$P_Beta,

    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. FIT ALL PHYLUM / TROPHIC SUBGROUPS
# ------------------------------------------------------------------------------

rows <- list()
qc_rows <- list()

for (group in GROUPS) {
  message("\nProcessing pTTPL subgroup fits: ", group)

  for (plot_id in PLOT_LIST) {

    weighted_file <- file.path(
      WEIGHTED_INPUT_ROOT, group, "UPARSE",
      paste0(plot_id, "_MPD_weighted_otu_All_abun.csv")
    )

    if (!file.exists(weighted_file)) {
      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        Group = group, Plot = plot_id,
        Weighted_file = weighted_file,
        File_exists = FALSE,
        Alignment_mode = NA_character_,
        N_input_taxa = NA_integer_,
        N_years_used = NA_integer_,
        stringsAsFactors = FALSE
      )
      next
    }

    mat0 <- read_weighted_matrix(weighted_file)
    ali <- align_matrix_to_metadata(
      mat0, plot_id, env,
      allow_positional = ALLOW_POSITIONAL_YEAR_FALLBACK
    )

    mat <- ali$mat
    meta <- ali$meta
    meta_one <- plot_metadata_one_row(meta)

    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      Group = group, Plot = plot_id,
      Weighted_file = weighted_file,
      File_exists = TRUE,
      Alignment_mode = ali$alignment_mode,
      N_input_taxa = nrow(mat),
      N_years_used = ncol(mat),
      stringsAsFactors = FALSE
    )

    # --------------------------------------------------------------------------
    # Phylum level
    # --------------------------------------------------------------------------
    if (group %in% PHYLUM_GROUPS) {
      taxonomy <- read_taxonomy(group, plot_id)

      if (!is.null(taxonomy) && "Phylum" %in% names(taxonomy)) {
        phy <- as.character(taxonomy$Phylum)
        names(phy) <- rownames(taxonomy)

        valid_label <- !is.na(phy) & nzchar(trimws(phy))
        phy <- phy[valid_label]

        tab <- sort(table(phy), decreasing = TRUE)
        candidate_phyla <- names(tab)[tab >= MIN_TAXA_PHYLUM]

        for (subgroup in candidate_phyla) {
          ids_membership <- names(phy)[phy == subgroup]
          ids_matched <- intersect(ids_membership, rownames(mat))

          if (length(ids_matched) == 0) next

          # IMPORTANT: no subgroup renormalization of MPD-weighted values.
          mat_sub <- mat[ids_matched, , drop = FALSE]

          fit_obj <- fit_subgroup(
            mat_sub, years = meta$year, min_taxa = MIN_TAXA_PHYLUM
          )

          rows[[length(rows) + 1]] <- make_result_row(
            group = group,
            level_type = "Phylum",
            subgroup = subgroup,
            meta_one = meta_one,
            n_membership = length(ids_membership),
            n_matched = length(ids_matched),
            fit_obj = fit_obj,
            min_taxa = MIN_TAXA_PHYLUM,
            alignment_mode = ali$alignment_mode
          )
        }
      }
    }

    # --------------------------------------------------------------------------
    # Trophic level
    # --------------------------------------------------------------------------
    if (group %in% names(TROPHIC_LEVELS)) {
      for (subgroup in TROPHIC_LEVELS[[group]]) {

        ids_membership <- read_trophic_membership(group, subgroup)
        ids_matched <- intersect(ids_membership, rownames(mat))

        if (length(ids_matched) == 0) next

        # Trophic files define membership only; pTTPL values are taken from the
        # whole-community MPD-weighted matrix without renormalization.
        mat_sub <- mat[ids_matched, , drop = FALSE]

        fit_obj <- fit_subgroup(
          mat_sub, years = meta$year, min_taxa = MIN_TAXA_TROPHIC
        )

        rows[[length(rows) + 1]] <- make_result_row(
          group = group,
          level_type = "Trophic",
          subgroup = subgroup,
          meta_one = meta_one,
          n_membership = length(ids_membership),
          n_matched = length(ids_matched),
          fit_obj = fit_obj,
          min_taxa = MIN_TAXA_TROPHIC,
          alignment_mode = ali$alignment_mode
        )
      }
    }
  }
}

results <- if (length(rows) > 0) {
  do.call(rbind, rows)
} else {
  data.frame()
}

qc <- if (length(qc_rows) > 0) {
  do.call(rbind, qc_rows)
} else {
  data.frame()
}

if (nrow(results) == 0) {
  stop("No subgroup pTTPL results were generated. Check input paths and OTU IDs.")
}

# ------------------------------------------------------------------------------
# 6. SAVE
# ------------------------------------------------------------------------------

out_csv <- file.path(
  OUTPUT_ROOT, "pTTPL_phylum_trophic_plot_level_full_precision.csv"
)
out_rds <- file.path(
  OUTPUT_ROOT, "pTTPL_phylum_trophic_plot_level_full_precision.rds"
)
qc_csv <- file.path(
  OUTPUT_ROOT, "pTTPL_phylum_trophic_input_QC.csv"
)

write_csv_full_precision(results, out_csv)
saveRDS(results, out_rds)
write_csv_full_precision(qc, qc_csv)

message("\nSaved:")
message("  ", out_csv)
message("  ", out_rds)
message("  ", qc_csv)
message("\nDone.")
