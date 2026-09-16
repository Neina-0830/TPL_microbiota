# ==============================================================================
# NWS pTTPL REVISION WORKFLOW — STEP 2 (categorical-primary aligned version)
# Complete beta_p hierarchical statistics from the full-precision pTTPL refit
#
# PRIMARY INPUT:
#   ./Taylor_revision_pTTPL/pTTPL_all_groups_plot_level_full_precision.rds
#
# PRIMARY MANUSCRIPT MODEL:
#   metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)
#
# Primary treatment factors:
#   WarmF   : Unwarmed / Warmed
#   PrecipF : Ambient / Half / Double
#   ClipF   : Unclipped / Clipped
#
# Numeric Warming / Precipitation / Clipping variables are retained only for
# the constrained continuous-precipitation sensitivity analysis and legacy
# compatibility.
#
# This script supplies:
#   - categorical hierarchical Type III tests
#   - coefficients, random-effect diagnostics, marginal/conditional R2
#   - block-fixed categorical mixed-model sensitivity
#   - ML comparison of categorical vs continuous-precipitation models
#   - 20 planned effect-specific contrasts, Holm within family
#   - four OVERALL factor effects from the SAME categorical LMM:
#       Overall Warming
#       Overall Half precipitation vs Ambient
#       Overall Double precipitation vs Ambient
#       Overall Clipping
#     with Holm correction across the four overall effects within Group x Metric
#   - 11 individual treatment-vs-Ctrl hierarchical contrasts for Table 1
#   - treatment-specific Cohen's d / Hedges' g for Table 1
#   - OVERALL Cohen's d / Hedges' g for Fig. 2a
#   - continuous precipitation slopes as a sensitivity analysis
#   - block-matched standardized sensitivity
#   - beta_p SE propagation with rma.mv, V = SE_Beta_p^2, using the same
#     categorical fixed-effect structure as the primary LMM
#   - beta_p/reference-variance warming-change co-variation
#   - optional beta_p association with Hellinger community temporal variability
#   - manuscript/reviewer-oriented compact tables and figures
#
# IMPORTANT:
#   Table 1 remains treatment-combination vs Ctrl.
#   Fig. 2a reports OVERALL factor effects.
#
# CR2 is intentionally NOT used in beta-SE propagation because Block is the
# highest-level dependency and there are only four top-level blocks.
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

INPUT_RDS <- "./results/04_pTTPL/pTTPL_all_groups_plot_level_full_precision.rds"
INPUT_CSV <- "./results/04_pTTPL/pTTPL_all_groups_plot_level_full_precision.csv"

OUTDIR <- "./results/05_pTTPL_factorial_statistics"
FIGDIR <- file.path(
  OUTDIR,
  "figures"
)

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIGDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Optional: use the already-computed taxonomic community temporal variability
# metric to evaluate whether beta_p is associated with the same whole-community
# compositional temporal variability measure.
RUN_CTV_ASSOCIATION <- TRUE

CTV_SOURCE_FILE <- "./results/01_tTTPL/UPARSE_all_groups_plot_level_revision_metrics.csv"

GROUP_ORDER <- c(
  "Archaea", "Bacteria", "Fungi", "Protists", "Nematoda"
)

TREATMENT_LEVELS <- c(
  "Ctrl", "W", "H", "D", "C",
  "WH", "WD", "WC", "HC", "DC",
  "WHC", "WDC"
)

# Primary reviewer metrics for pTTPL.
MAIN_METRICS <- c(
  "Beta_raw",
  "Beta_detrended",
  "LogVar_at_ref_raw",
  "LogVar_at_ref_detrended"
)

BETA_SE_PAIRS <- c(
  Beta_raw = "SE_Beta_raw",
  Beta_detrended = "SE_Beta_detrended"
)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "lme4",
  "lmerTest",
  "emmeans",
  "metafor"
)

missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install required packages first:\n  ",
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    "\n\nExample:\ninstall.packages(c(",
    paste(
      sprintf(
        '"%s"',
        missing_pkgs
      ),
      collapse = ", "
    ),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(metafor)
})

options(
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
)

emmeans::emm_options(
  lmer.df = "satterthwaite"
)

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

write_out <- function(
  x,
  filename
) {

  y <- as.data.frame(
    lapply(
      x,
      as_full_precision_character
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  write.table(
    y,
    file = file.path(
      OUTDIR,
      filename
    ),
    sep = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA",
    qmethod = "double"
  )
}

save_pdf_png <- function(
  p,
  stem,
  width,
  height
) {

  ggplot2::ggsave(
    file.path(
      FIGDIR,
      paste0(stem, ".pdf")
    ),
    p,
    width = width,
    height = height
  )

  ggplot2::ggsave(
    file.path(
      FIGDIR,
      paste0(stem, ".png")
    ),
    p,
    width = width,
    height = height,
    dpi = 300
  )
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    p <= 0.10 ~ "#",
    TRUE ~ "ns"
  )
}

# ------------------------------------------------------------------------------
# 3. PLANNED CONTRAST DEFINITIONS
# ------------------------------------------------------------------------------

contrast_def <- bind_rows(

  tibble(
    Family = "Warming",
    Contrast_order = 1:6,
    Target = c(
      "W", "WH", "WD",
      "WC", "WHC", "WDC"
    ),
    Reference = c(
      "Ctrl", "H", "D",
      "C", "HC", "DC"
    ),
    Context = c(
      "Ambient | Unclipped",
      "Half | Unclipped",
      "Double | Unclipped",
      "Ambient | Clipped",
      "Half | Clipped",
      "Double | Clipped"
    )
  ),

  tibble(
    Family = "Clipping",
    Contrast_order = 1:6,
    Target = c(
      "C", "WC", "HC",
      "DC", "WHC", "WDC"
    ),
    Reference = c(
      "Ctrl", "W", "H",
      "D", "WH", "WD"
    ),
    Context = c(
      "Unwarmed | Ambient",
      "Warmed | Ambient",
      "Unwarmed | Half",
      "Unwarmed | Double",
      "Warmed | Half",
      "Warmed | Double"
    )
  ),

  tibble(
    Family = "Half_precipitation",
    Contrast_order = 1:4,
    Target = c(
      "H", "WH", "HC", "WHC"
    ),
    Reference = c(
      "Ctrl", "W", "C", "WC"
    ),
    Context = c(
      "Unwarmed | Unclipped",
      "Warmed | Unclipped",
      "Unwarmed | Clipped",
      "Warmed | Clipped"
    )
  ),

  tibble(
    Family = "Double_precipitation",
    Contrast_order = 1:4,
    Target = c(
      "D", "WD", "DC", "WDC"
    ),
    Reference = c(
      "Ctrl", "W", "C", "WC"
    ),
    Context = c(
      "Unwarmed | Unclipped",
      "Warmed | Unclipped",
      "Unwarmed | Clipped",
      "Warmed | Clipped"
    )
  )
) %>%
  mutate(
    Contrast = paste0(
      Target,
      "-",
      Reference
    )
  )

# ------------------------------------------------------------------------------
# 4. READ EXACT pTTPL RESULTS
# ------------------------------------------------------------------------------

if (file.exists(INPUT_RDS)) {

  dat <- readRDS(
    INPUT_RDS
  )

  input_used <- INPUT_RDS

} else if (file.exists(INPUT_CSV)) {

  dat <- read.csv(
    INPUT_CSV,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  input_used <- INPUT_CSV

} else {

  stop(
    "Neither pTTPL input file exists:\n  ",
    INPUT_RDS,
    "\n  ",
    INPUT_CSV,
    "\nRun NWS_pTTPL_01_fit_full_precision.R first."
  )
}

required_core <- c(
  "Group", "Plot", "block",
  "WholePlotID", "Number",
  "Warming", "Precipitation",
  "Clipping", "Treatment",
  "Beta_raw", "SE_Beta_raw",
  "Beta_detrended",
  "SE_Beta_detrended",
  "LogVar_at_ref_raw",
  "LogVar_at_ref_detrended"
)

missing_core <- setdiff(
  required_core,
  names(dat)
)

if (length(missing_core) > 0) {
  stop(
    "pTTPL input is missing required columns: ",
    paste(
      missing_core,
      collapse = ", "
    )
  )
}

dat <- dat %>%
  mutate(
    Group = factor(
      as.character(Group),
      levels = GROUP_ORDER
    ),
    Plot = factor(
      as.character(Plot)
    ),
    Block = factor(
      as.character(block)
    ),
    Number = factor(
      as.character(Number)
    ),
    WholePlotID = factor(
      as.character(WholePlotID)
    ),
    BlockWholePlotID = interaction(
      Block,
      WholePlotID,
      drop = TRUE
    ),
    Treatment = factor(
      as.character(Treatment),
      levels = TREATMENT_LEVELS
    ),
    Warming = as.numeric(
      as.character(Warming)
    ),
    Precipitation = as.numeric(
      as.character(Precipitation)
    ),
    Clipping = as.numeric(
      as.character(Clipping)
    )
  )


# ------------------------------------------------------------------------------
# FACTOR VARIABLES FOR THE PRIMARY CATEGORICAL MODEL
# ------------------------------------------------------------------------------
# These factors define the primary experimental estimands. In particular,
# Half and Double precipitation are estimated independently rather than being
# constrained to a single linear precipitation slope.

dat$WarmF <- factor(
  ifelse(dat$Warming == 1, "Warmed", "Unwarmed"),
  levels = c("Unwarmed", "Warmed")
)

dat$PrecipF <- factor(
  dplyr::case_when(
    dat$Precipitation == -0.5 ~ "Half",
    dat$Precipitation == 0    ~ "Ambient",
    dat$Precipitation == 1    ~ "Double",
    TRUE ~ NA_character_
  ),
  levels = c("Ambient", "Half", "Double")
)

dat$ClipF <- factor(
  ifelse(dat$Clipping == 1, "Clipped", "Unclipped"),
  levels = c("Unclipped", "Clipped")
)

if (
  any(is.na(dat$WarmF)) ||
  any(is.na(dat$PrecipF)) ||
  any(is.na(dat$ClipF))
) {
  stop("Could not construct WarmF / PrecipF / ClipF from numeric treatment coding.")
}

# Check that the manuscript treatment label agrees with the factor combination.
dat$ExpectedTreatment <- dplyr::case_when(
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Ambient" & dat$ClipF == "Unclipped" ~ "Ctrl",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Ambient" & dat$ClipF == "Unclipped" ~ "W",
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Half"    & dat$ClipF == "Unclipped" ~ "H",
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Double"  & dat$ClipF == "Unclipped" ~ "D",
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Ambient" & dat$ClipF == "Clipped"   ~ "C",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Half"    & dat$ClipF == "Unclipped" ~ "WH",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Double"  & dat$ClipF == "Unclipped" ~ "WD",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Ambient" & dat$ClipF == "Clipped"   ~ "WC",
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Half"    & dat$ClipF == "Clipped"   ~ "HC",
  dat$WarmF == "Unwarmed" & dat$PrecipF == "Double"  & dat$ClipF == "Clipped"   ~ "DC",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Half"    & dat$ClipF == "Clipped"   ~ "WHC",
  dat$WarmF == "Warmed"   & dat$PrecipF == "Double"  & dat$ClipF == "Clipped"   ~ "WDC",
  TRUE ~ NA_character_
)

if (any(dat$ExpectedTreatment != as.character(dat$Treatment), na.rm = TRUE)) {
  bad_treatment <- dat %>%
    filter(ExpectedTreatment != as.character(Treatment)) %>%
    select(
      Group, Plot, Block, WholePlotID,
      WarmF, PrecipF, ClipF, Treatment,
      ExpectedTreatment
    )
  write_out(
    bad_treatment,
    "00_pTTPL_BAD_treatment_factor_mapping.csv"
  )
  stop(
    "Treatment labels do not agree with WarmF/PrecipF/ClipF. ",
    "See 00_pTTPL_BAD_treatment_factor_mapping.csv."
  )
}

dat$ExpectedTreatment <- NULL

if (any(is.na(dat$Group))) {
  stop("Unrecognized Group label in pTTPL input.")
}

if (any(is.na(dat$Treatment))) {
  stop("Unrecognized Treatment label in pTTPL input.")
}

if (!all(dat$Warming %in% c(0, 1))) {
  stop("Warming coding is not 0/1.")
}

if (!all(dat$Precipitation %in% c(-0.5, 0, 1))) {
  stop("Precipitation coding is not -0.5/0/1.")
}

if (!all(dat$Clipping %in% c(0, 1))) {
  stop("Clipping coding is not 0/1.")
}

MAIN_METRICS <- MAIN_METRICS[
  MAIN_METRICS %in% names(dat)
]

BETA_SE_PAIRS <- BETA_SE_PAIRS[
  names(BETA_SE_PAIRS) %in% names(dat) &
    unname(BETA_SE_PAIRS) %in% names(dat)
]

# ------------------------------------------------------------------------------
# 5. DESIGN QC
# ------------------------------------------------------------------------------

design_qc <- dat %>%
  group_by(Group) %>%
  summarise(
    N_rows = n(),
    N_plots = n_distinct(Plot),
    N_blocks = n_distinct(Block),
    N_wholeplots =
      n_distinct(BlockWholePlotID),
    N_treatments =
      n_distinct(Treatment),
    .groups = "drop"
  )

write_out(
  design_qc,
  "00_pTTPL_design_QC.csv"
)

nested_qc <- dat %>%
  distinct(
    Group,
    Block,
    BlockWholePlotID,
    Clipping
  ) %>%
  group_by(
    Group,
    BlockWholePlotID
  ) %>%
  summarise(
    N_blocks = n_distinct(Block),
    N_clipping_levels =
      n_distinct(Clipping),
    .groups = "drop"
  ) %>%
  mutate(
    Belongs_to_one_block =
      N_blocks == 1,
    Has_both_clipping_levels =
      N_clipping_levels == 2
  )

write_out(
  nested_qc,
  "01_pTTPL_nested_random_structure_QC.csv"
)

if (
  any(design_qc$N_rows != 48) ||
    any(design_qc$N_plots != 48) ||
    any(design_qc$N_blocks != 4) ||
    any(design_qc$N_wholeplots != 24) ||
    any(design_qc$N_treatments != 12)
) {
  warning(
    "At least one group does not match the expected 48-subplot, ",
    "24-whole-plot, 4-block, 12-treatment design. Inspect QC outputs."
  )
}

# ------------------------------------------------------------------------------
# 6. DESCRIPTIVE 12-TREATMENT SUMMARIES
# ------------------------------------------------------------------------------

descriptive_rows <- list()
idx <- 1

for (metric in MAIN_METRICS) {

  dd <- dat %>%
    filter(
      is.finite(.data[[metric]])
    ) %>%
    group_by(
      Group,
      Treatment
    ) %>%
    summarise(
      N = n(),
      Mean =
        mean(.data[[metric]]),
      SD =
        sd(.data[[metric]]),
      SE =
        SD / sqrt(N),
      .groups = "drop"
    ) %>%
    mutate(
      Metric = metric,
      .before = 2
    )

  descriptive_rows[[idx]] <- dd
  idx <- idx + 1
}

descriptive <- bind_rows(
  descriptive_rows
)

write_out(
  descriptive,
  "02_pTTPL_descriptive_12treatments.csv"
)

# ------------------------------------------------------------------------------
# 7. PRIMARY CATEGORICAL HIERARCHICAL MODEL
#    + CONTINUOUS-PRECIPITATION SENSITIVITY / MODEL COMPARISON
# ------------------------------------------------------------------------------

PRIMARY_MODEL_LABEL <-
  "metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)"

CONTINUOUS_MODEL_LABEL <-
  "metric ~ Warming * Precipitation * Clipping + (1 | Block/WholePlotID)"

BLOCK_FIXED_MODEL_LABEL <-
  "metric ~ WarmF * PrecipF * ClipF + Block + (1 | BlockWholePlotID)"

fit_primary_model <- function(
  d,
  metric
) {

  dd <- d %>%
    filter(
      is.finite(.data[[metric]])
    ) %>%
    droplevels()

  if (nrow(dd) < 20) {
    return(list(
      fit = NULL,
      continuous = NULL,
      fit_ml = NULL,
      continuous_ml = NULL,
      data = dd,
      error = "Too few observations."
    ))
  }

  form_cat <- as.formula(
    paste0(
      metric,
      " ~ WarmF * PrecipF * ClipF + ",
      "(1 | Block/WholePlotID)"
    )
  )

  form_num <- as.formula(
    paste0(
      metric,
      " ~ Warming * Precipitation * Clipping + ",
      "(1 | Block/WholePlotID)"
    )
  )

  # REML categorical model = PRIMARY inference.
  fit_cat <- tryCatch(
    lmerTest::lmer(
      form_cat,
      data = dd,
      REML = TRUE
    ),
    error = function(e) e
  )

  # REML continuous-gradient model = sensitivity only.
  fit_num <- tryCatch(
    lmerTest::lmer(
      form_num,
      data = dd,
      REML = TRUE
    ),
    error = function(e) e
  )

  # ML fits are used ONLY for comparing fixed-effect structures.
  fit_cat_ml <- tryCatch(
    lme4::lmer(
      form_cat,
      data = dd,
      REML = FALSE
    ),
    error = function(e) e
  )

  fit_num_ml <- tryCatch(
    lme4::lmer(
      form_num,
      data = dd,
      REML = FALSE
    ),
    error = function(e) e
  )

  errors <- character(0)

  if (inherits(fit_cat, "error")) {
    errors <- c(
      errors,
      paste0(
        "categorical REML: ",
        conditionMessage(fit_cat)
      )
    )
    fit_cat <- NULL
  }

  if (inherits(fit_num, "error")) {
    errors <- c(
      errors,
      paste0(
        "continuous REML: ",
        conditionMessage(fit_num)
      )
    )
    fit_num <- NULL
  }

  if (inherits(fit_cat_ml, "error")) {
    errors <- c(
      errors,
      paste0(
        "categorical ML: ",
        conditionMessage(fit_cat_ml)
      )
    )
    fit_cat_ml <- NULL
  }

  if (inherits(fit_num_ml, "error")) {
    errors <- c(
      errors,
      paste0(
        "continuous ML: ",
        conditionMessage(fit_num_ml)
      )
    )
    fit_num_ml <- NULL
  }

  list(
    fit = fit_cat,
    continuous = fit_num,
    fit_ml = fit_cat_ml,
    continuous_ml = fit_num_ml,
    data = dd,
    error =
      if (length(errors) == 0) {
        NA_character_
      } else {
        paste(
          errors,
          collapse = " | "
        )
      }
  )
}


fit_blockfixed_model <- function(
  d,
  metric
) {

  form <- as.formula(
    paste0(
      metric,
      " ~ WarmF * PrecipF * ClipF + ",
      "Block + (1 | BlockWholePlotID)"
    )
  )

  lmerTest::lmer(
    form,
    data = d,
    REML = TRUE
  )
}


# Nakagawa-style R2 for the present Gaussian random-intercept LMM.
nakagawa_r2_gaussian <- function(model) {

  if (is.null(model)) {
    return(c(
      R2_marginal = NA_real_,
      R2_conditional = NA_real_
    ))
  }

  fixed_pred <- tryCatch(
    as.numeric(
      stats::predict(
        model,
        re.form = NA
      )
    ),
    error = function(e)
      rep(
        NA_real_,
        stats::nobs(model)
      )
  )

  fixed_var <-
    if (
      sum(
        is.finite(fixed_pred)
      ) > 1
    ) {
      stats::var(
        fixed_pred,
        na.rm = TRUE
      )
    } else {
      NA_real_
    }

  vc <- as.data.frame(
    lme4::VarCorr(model)
  )

  random_var <- sum(
    vc$vcov[
      vc$grp != "Residual"
    ],
    na.rm = TRUE
  )

  residual_var <- vc %>%
    filter(
      grp == "Residual"
    ) %>%
    pull(vcov)

  if (length(residual_var) == 0) {
    residual_var <- NA_real_
  }

  residual_var <- residual_var[1]

  total_var <-
    fixed_var +
    random_var +
    residual_var

  c(
    R2_marginal =
      ifelse(
        is.finite(total_var) &&
          total_var > 0,
        fixed_var / total_var,
        NA_real_
      ),
    R2_conditional =
      ifelse(
        is.finite(total_var) &&
          total_var > 0,
        (fixed_var + random_var) /
          total_var,
        NA_real_
      )
  )
}


extract_diag <- function(
  model,
  group,
  metric,
  data_used
) {

  if (is.null(model)) {
    return(tibble(
      Group = group,
      Metric = metric,
      Model = PRIMARY_MODEL_LABEL,
      N = nrow(data_used),
      N_blocks =
        n_distinct(data_used$Block),
      N_nested_wholeplots =
        n_distinct(
          data_used$BlockWholePlotID
        ),
      Singular = NA,
      Block_variance = NA_real_,
      WholePlot_within_Block_variance =
        NA_real_,
      Residual_variance = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      Residual_mean = NA_real_,
      Residual_SD = NA_real_,
      Max_abs_standardized_residual =
        NA_real_,
      N_abs_standardized_residual_gt3 =
        NA_integer_,
      QQ_correlation = NA_real_,
      AbsResidual_vs_Fitted_P =
        NA_real_
    ))
  }

  vc <- as.data.frame(
    VarCorr(model)
  )

  block_var <- vc %>%
    filter(
      grp == "Block"
    ) %>%
    pull(vcov)

  if (length(block_var) == 0) {
    block_var <- NA_real_
  }

  wp_var <- vc %>%
    filter(
      grepl(
        "Block",
        grp
      ) &
        grepl(
          "WholePlotID",
          grp
        )
    ) %>%
    pull(vcov)

  if (length(wp_var) == 0) {
    wp_var <- NA_real_
  }

  resid_var <- vc %>%
    filter(
      grp == "Residual"
    ) %>%
    pull(vcov)

  if (length(resid_var) == 0) {
    resid_var <- NA_real_
  }

  r2 <- nakagawa_r2_gaussian(
    model
  )

  residual <- stats::residuals(
    model
  )

  fitted_value <- stats::fitted(
    model
  )

  std_residual <-
    residual /
    stats::sigma(model)

  qq_theory <-
    stats::qnorm(
      stats::ppoints(
        length(std_residual)
      )
    )

  qq_observed <-
    sort(std_residual)

  qq_cor <-
    if (
      length(std_residual) > 2 &&
      is.finite(
        stats::sd(
          std_residual
        )
      ) &&
      stats::sd(
        std_residual
      ) > 0
    ) {
      stats::cor(
        qq_theory,
        qq_observed
      )
    } else {
      NA_real_
    }

  hetero_p <- tryCatch(
    summary(
      stats::lm(
        abs(std_residual) ~
          fitted_value
      )
    )$coefficients[2, 4],
    error = function(e)
      NA_real_
  )

  tibble(
    Group = group,
    Metric = metric,
    Model = PRIMARY_MODEL_LABEL,
    N = nrow(data_used),
    N_blocks =
      n_distinct(data_used$Block),
    N_nested_wholeplots =
      n_distinct(
        data_used$BlockWholePlotID
      ),
    Singular =
      lme4::isSingular(
        model,
        tol = 1e-6
      ),
    Block_variance =
      block_var[1],
    WholePlot_within_Block_variance =
      wp_var[1],
    Residual_variance =
      resid_var[1],
    R2_marginal =
      unname(
        r2["R2_marginal"]
      ),
    R2_conditional =
      unname(
        r2["R2_conditional"]
      ),
    Residual_mean =
      mean(
        residual,
        na.rm = TRUE
      ),
    Residual_SD =
      sd(
        residual,
        na.rm = TRUE
      ),
    Max_abs_standardized_residual =
      max(
        abs(std_residual),
        na.rm = TRUE
      ),
    N_abs_standardized_residual_gt3 =
      sum(
        abs(std_residual) > 3,
        na.rm = TRUE
      ),
    QQ_correlation =
      qq_cor,
    AbsResidual_vs_Fitted_P =
      hetero_p
  )
}


compare_cat_vs_continuous <- function(
  m_cat_ml,
  m_num_ml,
  group,
  metric
) {

  if (
    is.null(m_cat_ml) ||
    is.null(m_num_ml)
  ) {
    return(tibble(
      Group = group,
      Metric = metric,
      AIC_categorical = NA_real_,
      AIC_continuous = NA_real_,
      Delta_AIC_cat_minus_cont = NA_real_,
      BIC_categorical = NA_real_,
      BIC_continuous = NA_real_,
      Delta_BIC_cat_minus_cont = NA_real_,
      LogLik_categorical = NA_real_,
      LogLik_continuous = NA_real_,
      LR_Chisq = NA_real_,
      LR_df = NA_real_,
      LR_P = NA_real_
    ))
  }

  aic_cat <- AIC(
    m_cat_ml
  )
  aic_num <- AIC(
    m_num_ml
  )
  bic_cat <- BIC(
    m_cat_ml
  )
  bic_num <- BIC(
    m_num_ml
  )

  ll_cat <- as.numeric(
    logLik(m_cat_ml)
  )
  ll_num <- as.numeric(
    logLik(m_num_ml)
  )

  lr <- tryCatch(
    stats::anova(
      m_num_ml,
      m_cat_ml,
      refit = FALSE
    ),
    error = function(e) NULL
  )

  lr_chisq <- lr_df <-
    lr_p <- NA_real_

  if (
    !is.null(lr) &&
    nrow(lr) >= 2
  ) {
    if ("Chisq" %in% names(lr)) {
      lr_chisq <-
        as.numeric(
          lr$Chisq[2]
        )
    }
    if ("Chi Df" %in% names(lr)) {
      lr_df <-
        as.numeric(
          lr[["Chi Df"]][2]
        )
    }
    if ("Pr(>Chisq)" %in% names(lr)) {
      lr_p <-
        as.numeric(
          lr[["Pr(>Chisq)"]][2]
        )
    }
  }

  tibble(
    Group = group,
    Metric = metric,
    AIC_categorical =
      aic_cat,
    AIC_continuous =
      aic_num,
    Delta_AIC_cat_minus_cont =
      aic_cat - aic_num,
    BIC_categorical =
      bic_cat,
    BIC_continuous =
      bic_num,
    Delta_BIC_cat_minus_cont =
      bic_cat - bic_num,
    LogLik_categorical =
      ll_cat,
    LogLik_continuous =
      ll_num,
    LR_Chisq =
      lr_chisq,
    LR_df =
      lr_df,
    LR_P =
      lr_p
  )
}


primary_models <- list()
type3_rows <- list()
coef_rows <- list()
diag_rows <- list()
fixed_rows <- list()
model_comparison_rows <- list()

idx <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    ) %>%
    droplevels()

  for (metric in MAIN_METRICS) {

    ans <- fit_primary_model(
      dg,
      metric
    )

    key <- paste(
      g,
      metric,
      sep = "::"
    )

    primary_models[[key]] <- ans

    if (!is.null(ans$fit)) {

      m <- ans$fit

      an <- as.data.frame(
        stats::anova(
          m,
          type = "III",
          ddf =
            "Satterthwaite"
        )
      )

      an$Term <- rownames(an)
      rownames(an) <- NULL

      type3_rows[[idx]] <- an %>%
        mutate(
          Group = g,
          Metric = metric,
          Model =
            PRIMARY_MODEL_LABEL,
          .before = 1
        )

      co <- as.data.frame(
        summary(m)$coefficients
      )

      co$Term <- rownames(co)
      rownames(co) <- NULL

      coef_rows[[idx]] <- co %>%
        mutate(
          Group = g,
          Metric = metric,
          Model =
            PRIMARY_MODEL_LABEL,
          .before = 1
        )

      diag_rows[[idx]] <-
        extract_diag(
          m,
          g,
          metric,
          ans$data
        )

      mf <- tryCatch(
        fit_blockfixed_model(
          ans$data,
          metric
        ),
        error = function(e) e
      )

      if (!inherits(mf, "error")) {

        af <- as.data.frame(
          stats::anova(
            mf,
            type = "III",
            ddf =
              "Satterthwaite"
          )
        )

        af$Term <- rownames(af)
        rownames(af) <- NULL

        fixed_rows[[idx]] <- af %>%
          mutate(
            Group = g,
            Metric = metric,
            Model =
              BLOCK_FIXED_MODEL_LABEL,
            Singular =
              lme4::isSingular(
                mf,
                tol = 1e-6
              ),
            .before = 1
          )
      }

    } else {

      diag_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        Model =
          PRIMARY_MODEL_LABEL,
        N = nrow(ans$data),
        N_blocks =
          n_distinct(ans$data$Block),
        N_nested_wholeplots =
          n_distinct(
            ans$data$BlockWholePlotID
          ),
        Singular = NA,
        Block_variance =
          NA_real_,
        WholePlot_within_Block_variance =
          NA_real_,
        Residual_variance =
          NA_real_,
        R2_marginal =
          NA_real_,
        R2_conditional =
          NA_real_,
        Residual_mean =
          NA_real_,
        Residual_SD =
          NA_real_,
        Max_abs_standardized_residual =
          NA_real_,
        N_abs_standardized_residual_gt3 =
          NA_integer_,
        QQ_correlation =
          NA_real_,
        AbsResidual_vs_Fitted_P =
          NA_real_,
        Error = ans$error
      )
    }

    model_comparison_rows[[idx]] <-
      compare_cat_vs_continuous(
        ans$fit_ml,
        ans$continuous_ml,
        g,
        metric
      )

    idx <- idx + 1
  }
}


type3_out <- bind_rows(
  type3_rows
)

coef_out <- bind_rows(
  coef_rows
)

diag_out <- bind_rows(
  diag_rows
)

fixed_out <- bind_rows(
  fixed_rows
)

model_comparison_out <- bind_rows(
  model_comparison_rows
)


write_out(
  type3_out,
  "03_pTTPL_PRIMARY_categorical_hierarchical_TypeIII.csv"
)

write_out(
  coef_out,
  "03_pTTPL_PRIMARY_categorical_model_coefficients.csv"
)

write_out(
  diag_out,
  "03_pTTPL_PRIMARY_categorical_model_diagnostics.csv"
)

write_out(
  fixed_out,
  "03_pTTPL_BLOCK_FIXED_categorical_sensitivity_TypeIII.csv"
)

write_out(
  model_comparison_out,
  "03_pTTPL_MODEL_COMPARISON_categorical_vs_continuous.csv"
)


# Backward-compatible filenames. These now contain the categorical-primary
# results; the filename is retained only so downstream manuscript scripts do not
# fail when re-run.
write_out(
  type3_out,
  "03_pTTPL_PRIMARY_hierarchical_factorial_TypeIII.csv"
)

write_out(
  coef_out,
  "03_pTTPL_PRIMARY_standard_model_coefficients.csv"
)

write_out(
  diag_out,
  "03_pTTPL_PRIMARY_hierarchical_model_diagnostics.csv"
)

write_out(
  fixed_out,
  "03_pTTPL_blockfixed_sensitivity_TypeIII.csv"
)

# ------------------------------------------------------------------------------
# 8. MODEL-BASED PLANNED CONTRASTS FROM THE CATEGORICAL PRIMARY MODEL
# ------------------------------------------------------------------------------

make_treatment <- function(
  warm,
  precip,
  clip
) {

  key <- paste(
    as.character(warm),
    as.character(precip),
    as.character(clip),
    sep = "|"
  )

  map <- c(
    "Unwarmed|Ambient|Unclipped" = "Ctrl",
    "Warmed|Ambient|Unclipped"   = "W",
    "Unwarmed|Half|Unclipped"    = "H",
    "Unwarmed|Double|Unclipped"  = "D",
    "Unwarmed|Ambient|Clipped"   = "C",
    "Warmed|Half|Unclipped"      = "WH",
    "Warmed|Double|Unclipped"    = "WD",
    "Warmed|Ambient|Clipped"     = "WC",
    "Unwarmed|Half|Clipped"      = "HC",
    "Unwarmed|Double|Clipped"    = "DC",
    "Warmed|Half|Clipped"        = "WHC",
    "Warmed|Double|Clipped"      = "WDC"
  )

  unname(
    map[key]
  )
}


make_emm_weight <- function(
  emm_levels,
  target,
  reference
) {

  w <- rep(
    0,
    length(emm_levels)
  )

  names(w) <- emm_levels

  if (!all(
    c(
      target,
      reference
    ) %in% emm_levels
  )) {
    stop(
      "Missing EMM treatment level(s): ",
      paste(
        setdiff(
          c(target, reference),
          emm_levels
        ),
        collapse = ", "
      )
    )
  }

  w[target] <- 1
  w[reference] <- -1

  w
}


make_average_weight <- function(x) {
  Reduce(
    "+",
    x
  ) /
    length(x)
}


extract_primary_contrasts <- function(
  model,
  group,
  metric
) {

  # Full 2 x 3 x 2 categorical EMM grid.
  emm <- emmeans::emmeans(
    model,
    specs =
      ~ WarmF *
      PrecipF *
      ClipF,
    weights = "equal",
    lmer.df =
      "satterthwaite"
  )

  emm_grid <- as.data.frame(
    emm
  )

  emm_levels <- make_treatment(
    emm_grid$WarmF,
    emm_grid$PrecipF,
    emm_grid$ClipF
  )

  if (
    any(is.na(emm_levels)) ||
      anyDuplicated(emm_levels) > 0 ||
      length(emm_levels) != 12
  ) {
    stop(
      "Could not reconstruct the unique 12-treatment categorical EMM grid."
    )
  }

  individual_rows <- list()
  weights_by_family <- list()

  for (
    fam in unique(
      contrast_def$Family
    )
  ) {

    defs <- contrast_def %>%
      filter(
        Family == fam
      ) %>%
      arrange(
        Contrast_order
      )

    fam_weights <- list()

    for (
      j in seq_len(
        nrow(defs)
      )
    ) {

      w <- make_emm_weight(
        emm_levels,
        defs$Target[j],
        defs$Reference[j]
      )

      fam_weights[[j]] <- w

      cc <- emmeans::contrast(
        emm,
        method = setNames(
          list(w),
          defs$Contrast[j]
        ),
        adjust = "none"
      )

      ss <- as.data.frame(
        summary(
          cc,
          infer = c(
            TRUE,
            TRUE
          ),
          adjust = "none"
        )
      )

      individual_rows[[
        length(
          individual_rows
        ) + 1
      ]] <- tibble(
        Group = group,
        Metric = metric,
        Family = fam,
        Contrast_order =
          defs$Contrast_order[j],
        Contrast =
          defs$Contrast[j],
        Target =
          defs$Target[j],
        Reference =
          defs$Reference[j],
        Context =
          defs$Context[j],
        Model_basis =
          PRIMARY_MODEL_LABEL,
        Estimate =
          ss$estimate[1],
        SE =
          ss$SE[1],
        df =
          ss$df[1],
        CI_low =
          ss$lower.CL[1],
        CI_high =
          ss$upper.CL[1],
        t_ratio =
          ss$t.ratio[1],
        P_raw =
          ss$p.value[1]
      )
    }

    weights_by_family[[fam]] <-
      fam_weights
  }

  individual <- bind_rows(
    individual_rows
  ) %>%
    group_by(
      Group,
      Metric,
      Family
    ) %>%
    mutate(
      P_Holm_within_family =
        p.adjust(
          P_raw,
          method = "holm"
        ),
      Significance_Holm =
        sig_label(
          P_Holm_within_family
        )
    ) %>%
    ungroup() %>%
    group_by(
      Group,
      Metric
    ) %>%
    mutate(
      P_Holm_all20 =
        p.adjust(
          P_raw,
          method = "holm"
        )
    ) %>%
    ungroup()


  # ---------------------------------------------------------------------------
  # OVERALL factor effects
  #
  # Warming:
  #   mean of six simple warming effects across precipitation x clipping.
  #
  # Half:
  #   mean of four Half - Ambient effects across warming x clipping.
  #
  # Double:
  #   mean of four Double - Ambient effects across warming x clipping.
  #
  # Clipping:
  #   mean of six clipping effects across warming x precipitation.
  #
  # The four raw P values are Holm-adjusted together within Group x Metric.
  # ---------------------------------------------------------------------------

  overall_rows <- list()

  for (
    fam in names(
      weights_by_family
    )
  ) {

    wavg <- make_average_weight(
      weights_by_family[[fam]]
    )

    label <- paste0(
      "Overall_",
      fam
    )

    cc <- emmeans::contrast(
      emm,
      method = setNames(
        list(wavg),
        label
      ),
      adjust = "none"
    )

    ss <- as.data.frame(
      summary(
        cc,
        infer = c(
          TRUE,
          TRUE
        ),
        adjust = "none"
      )
    )

    overall_rows[[
      length(overall_rows) + 1
    ]] <- tibble(
      Group = group,
      Metric = metric,
      Family = fam,
      Contrast = label,
      Model_basis =
        PRIMARY_MODEL_LABEL,
      Estimate =
        ss$estimate[1],
      SE =
        ss$SE[1],
      df =
        ss$df[1],
      CI_low =
        ss$lower.CL[1],
      CI_high =
        ss$upper.CL[1],
      t_ratio =
        ss$t.ratio[1],
      P_raw =
        ss$p.value[1]
    )
  }

  overall <- bind_rows(
    overall_rows
  ) %>%
    group_by(
      Group,
      Metric
    ) %>%
    mutate(
      P_Holm_four_overall =
        p.adjust(
          P_raw,
          method = "holm"
        ),
      Significance_raw =
        sig_label(
          P_raw
        ),
      Significance_Holm =
        sig_label(
          P_Holm_four_overall
        )
    ) %>%
    ungroup()

  list(
    individual = individual,
    overall = overall,
    emm = emm,
    emm_levels = emm_levels
  )
}


primary_contrast_rows <- list()
primary_overall_rows <- list()

idx <- 1

for (key in names(primary_models)) {

  ans <- primary_models[[key]]

  if (is.null(ans$fit)) {
    next
  }

  pieces <- strsplit(
    key,
    "::",
    fixed = TRUE
  )[[1]]

  g <- pieces[1]
  metric <- pieces[2]

  rr <- tryCatch(
    extract_primary_contrasts(
      ans$fit,
      g,
      metric
    ),
    error = function(e) e
  )

  if (!inherits(rr, "error")) {

    primary_contrast_rows[[idx]] <-
      rr$individual

    primary_overall_rows[[idx]] <-
      rr$overall

    idx <- idx + 1

  } else {

    warning(
      "Categorical contrast extraction failed for ",
      key,
      ": ",
      conditionMessage(rr)
    )
  }
}


primary_contrasts <- bind_rows(
  primary_contrast_rows
)

primary_overall <- bind_rows(
  primary_overall_rows
)


write_out(
  primary_contrasts,
  "04_pTTPL_PRIMARY_categorical_effect_specific_contrasts_Holm.csv"
)

write_out(
  primary_overall,
  "05_pTTPL_PRIMARY_categorical_overall_effects_Holm.csv"
)

# Backward-compatible filenames for manuscript/reviewer code that already reads
# these names. Their contents now come from the categorical primary model.
write_out(
  primary_contrasts,
  "04_pTTPL_PRIMARY_effect_specific_contrasts_Holm.csv"
)

write_out(
  primary_overall,
  "05_pTTPL_PRIMARY_overall_directional_effects.csv"
)

# ------------------------------------------------------------------------------
# 9. CONTINUOUS PRECIPITATION-GRADIENT SENSITIVITY
# ------------------------------------------------------------------------------
# This section is NOT the primary precipitation analysis. It summarizes the
# constrained numeric precipitation model (-0.5, 0, +1) for sensitivity and
# comparison with the categorical primary model.
# ------------------------------------------------------------------------------

precip_rows <- list()
idx <- 1

for (key in names(primary_models)) {

  ans <- primary_models[[key]]

  if (is.null(ans$continuous)) {
    next
  }

  pieces <- strsplit(
    key,
    "::",
    fixed = TRUE
  )[[1]]

  g <- pieces[1]
  metric <- pieces[2]
  model <- ans$continuous

  tr <- tryCatch(
    emmeans::emtrends(
      model,
      specs =
        ~ Warming * Clipping,
      var = "Precipitation",
      at = list(
        Warming = c(0, 1),
        Clipping = c(0, 1)
      ),
      lmer.df =
        "satterthwaite"
    ),
    error = function(e) e
  )

  if (inherits(tr, "error")) {
    warning(
      "Continuous precipitation emtrends failed for ",
      key,
      ": ",
      conditionMessage(tr)
    )
    next
  }

  ss <- as.data.frame(
    summary(
      tr,
      infer = c(
        TRUE,
        TRUE
      ),
      adjust = "none"
    )
  )

  trend_col <- grep(
    "\\.trend$",
    names(ss),
    value = TRUE
  )

  if (length(trend_col) != 1) {
    stop(
      "Could not identify precipitation trend column for ",
      key
    )
  }

  context <- tibble(
    Group = g,
    Metric = metric,
    Scope =
      "Context_specific",
    Warming =
      ss$Warming,
    Clipping =
      ss$Clipping,
    Context = paste0(
      ifelse(
        ss$Warming == 1,
        "Warmed",
        "Unwarmed"
      ),
      " | ",
      ifelse(
        ss$Clipping == 1,
        "Clipped",
        "Unclipped"
      )
    ),
    Precipitation_slope =
      ss[[trend_col]],
    SE = ss$SE,
    df = ss$df,
    CI_low = ss$lower.CL,
    CI_high = ss$upper.CL,
    t_ratio = ss$t.ratio,
    P_raw = ss$p.value
  )

  w <- rep(
    1 / nrow(ss),
    nrow(ss)
  )

  ov <- emmeans::contrast(
    tr,
    method = list(
      Overall_precipitation_slope =
        w
    ),
    adjust = "none"
  )

  oo <- as.data.frame(
    summary(
      ov,
      infer = c(
        TRUE,
        TRUE
      ),
      adjust = "none"
    )
  )

  overall <- tibble(
    Group = g,
    Metric = metric,
    Scope = "Overall",
    Warming = NA_real_,
    Clipping = NA_real_,
    Context =
      "Equal-weight average across warming x clipping contexts",
    Precipitation_slope =
      oo$estimate[1],
    SE =
      oo$SE[1],
    df =
      oo$df[1],
    CI_low =
      oo$lower.CL[1],
    CI_high =
      oo$upper.CL[1],
    t_ratio =
      oo$t.ratio[1],
    P_raw =
      oo$p.value[1]
  )

  precip_rows[[idx]] <- bind_rows(
    context,
    overall
  )

  idx <- idx + 1
}


precip_slopes <- bind_rows(
  precip_rows
) %>%
  group_by(
    Group,
    Metric,
    Scope
  ) %>%
  mutate(
    P_Holm_within_scope =
      p.adjust(
        P_raw,
        method = "holm"
      )
  ) %>%
  ungroup()


write_out(
  precip_slopes,
  "05_pTTPL_SENSITIVITY_continuous_precipitation_slopes.csv"
)

# Backward-compatible filename.
write_out(
  precip_slopes,
  "05_pTTPL_PRIMARY_continuous_precipitation_slopes.csv"
)

# ------------------------------------------------------------------------------
# 10. ALL 11 TREATMENTS VS CTRL FROM THE SAME HIERARCHICAL MODEL
# ------------------------------------------------------------------------------

vsctrl_rows <- list()
idx <- 1

for (key in names(primary_models)) {

  ans <- primary_models[[key]]

  if (is.null(ans$fit)) {
    next
  }

  pieces <- strsplit(
    key,
    "::",
    fixed = TRUE
  )[[1]]

  g <- pieces[1]
  metric <- pieces[2]

  rr <- extract_primary_contrasts(
    ans$fit,
    g,
    metric
  )

  emm <- rr$emm
  emm_levels <- rr$emm_levels

  for (
    trt in setdiff(
      TREATMENT_LEVELS,
      "Ctrl"
    )
  ) {

    w <- make_emm_weight(
      emm_levels,
      trt,
      "Ctrl"
    )

    label <- paste0(
      trt,
      "-Ctrl"
    )

    cc <- emmeans::contrast(
      emm,
      method = setNames(
        list(w),
        label
      ),
      adjust = "none"
    )

    ss <- as.data.frame(
      summary(
        cc,
        infer = c(
          TRUE,
          TRUE
        ),
        adjust = "none"
      )
    )

    vsctrl_rows[[idx]] <- tibble(
      Group = g,
      Metric = metric,
      Treatment = trt,
      Reference = "Ctrl",
      Contrast = label,
      Estimate =
        ss$estimate[1],
      SE =
        ss$SE[1],
      df =
        ss$df[1],
      CI_low =
        ss$lower.CL[1],
      CI_high =
        ss$upper.CL[1],
      t_ratio =
        ss$t.ratio[1],
      P_raw =
        ss$p.value[1]
    )

    idx <- idx + 1
  }
}

vsctrl <- bind_rows(
  vsctrl_rows
) %>%
  group_by(
    Group,
    Metric
  ) %>%
  mutate(
    P_Holm_11_vsCtrl =
      p.adjust(
        P_raw,
        method = "holm"
      ),
    Significance_Holm =
      sig_label(
        P_Holm_11_vsCtrl
      )
  ) %>%
  ungroup()

write_out(
  vsctrl,
  "06_pTTPL_PRIMARY_all_treatments_vs_Ctrl.csv"
)

# ------------------------------------------------------------------------------
# 11. DESCRIPTIVE COHEN'S d / HEDGES' g
# ------------------------------------------------------------------------------

standardized_difference <- function(
  x_target,
  x_reference
) {

  x_target <- x_target[
    is.finite(x_target)
  ]

  x_reference <- x_reference[
    is.finite(x_reference)
  ]

  nt <- length(x_target)
  nr <- length(x_reference)

  if (nt < 2 || nr < 2) {
    return(
      tibble(
        N_target = nt,
        N_reference = nr,
        Mean_target = NA_real_,
        Mean_reference = NA_real_,
        Raw_difference = NA_real_,
        Cohen_d_pooled = NA_real_,
        Hedges_g = NA_real_
      )
    )
  }

  st <- sd(x_target)
  sr <- sd(x_reference)

  pooled_sd <- sqrt(
    (
      (nt - 1) * st^2 +
        (nr - 1) * sr^2
    ) /
      (nt + nr - 2)
  )

  d <- if (
    is.finite(pooled_sd) &&
      pooled_sd > 0
  ) {
    (
      mean(x_target) -
        mean(x_reference)
    ) /
      pooled_sd
  } else {
    NA_real_
  }

  J <- 1 - 3 /
    (
      4 * (nt + nr) -
        9
    )

  tibble(
    N_target = nt,
    N_reference = nr,
    Mean_target =
      mean(x_target),
    Mean_reference =
      mean(x_reference),
    Raw_difference =
      mean(x_target) -
      mean(x_reference),
    Cohen_d_pooled = d,
    Hedges_g = J * d
  )
}

cohen_rows <- list()
idx <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    )

  for (metric in MAIN_METRICS) {

    for (
      j in seq_len(
        nrow(contrast_def)
      )
    ) {

      cd <- contrast_def[j, ]

      xt <- dg %>%
        filter(
          Treatment ==
            cd$Target
        ) %>%
        pull(
          all_of(metric)
        )

      xr <- dg %>%
        filter(
          Treatment ==
            cd$Reference
        ) %>%
        pull(
          all_of(metric)
        )

      cohen_rows[[idx]] <- bind_cols(
        tibble(
          Group = g,
          Metric = metric,
          Family =
            cd$Family,
          Contrast_order =
            cd$Contrast_order,
          Contrast =
            cd$Contrast,
          Target =
            cd$Target,
          Reference =
            cd$Reference,
          Context =
            cd$Context
        ),
        standardized_difference(
          xt,
          xr
        )
      )

      idx <- idx + 1
    }
  }
}

cohen_out <- bind_rows(
  cohen_rows
)

write_out(
  cohen_out,
  "07_pTTPL_Cohen_d_Hedges_g_effect_specific.csv"
)

# Descriptive values for the manuscript's treatment-vs-Ctrl display.
cohen_vsctrl_rows <- list()
idx <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    )

  for (metric in MAIN_METRICS) {

    xr <- dg %>%
      filter(
        Treatment == "Ctrl"
      ) %>%
      pull(
        all_of(metric)
      )

    for (
      trt in setdiff(
        TREATMENT_LEVELS,
        "Ctrl"
      )
    ) {

      xt <- dg %>%
        filter(
          Treatment == trt
        ) %>%
        pull(
          all_of(metric)
        )

      cohen_vsctrl_rows[[idx]] <- bind_cols(
        tibble(
          Group = g,
          Metric = metric,
          Treatment = trt,
          Reference = "Ctrl",
          Contrast = paste0(
            trt,
            "-Ctrl"
          )
        ),
        standardized_difference(
          xt,
          xr
        )
      )

      idx <- idx + 1
    }
  }
}

cohen_vsctrl <- bind_rows(
  cohen_vsctrl_rows
)

write_out(
  cohen_vsctrl,
  "07_pTTPL_Cohen_d_Hedges_g_vs_Ctrl.csv"
)


# ------------------------------------------------------------------------------
# 11B. OVERALL COHEN'S d / HEDGES' g FOR FIG. 2a
# ------------------------------------------------------------------------------
# These standardized effect sizes correspond to the SAME factor-level questions
# as the categorical-LMM overall contrasts:
#
# Warming: Warmed vs Unwarmed, pooling equally represented precipitation x
#          clipping contexts (24 vs 24 observations in the balanced design).
# Half:    Half vs Ambient, pooling warming x clipping contexts (16 vs 16).
# Double:  Double vs Ambient, pooling warming x clipping contexts (16 vs 16).
# Clipping: Clipped vs Unclipped, pooling warming x precipitation contexts
#           (24 vs 24).
#
# These Cohen/Hedges values are DESCRIPTIVE only. Formal inference is supplied by
# the categorical hierarchical LMM overall contrasts.

overall_cohen_rows <- list()
idx <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    ) %>%
    droplevels()

  for (metric in MAIN_METRICS) {

    if (!metric %in% names(dg)) {
      next
    }

    definitions <- list(
      list(
        Family = "Warming",
        Effect = "Overall_Warming",
        Target = "Warmed",
        Reference = "Unwarmed",
        target_values =
          dg %>%
          filter(
            WarmF == "Warmed"
          ) %>%
          pull(
            all_of(metric)
          ),
        reference_values =
          dg %>%
          filter(
            WarmF == "Unwarmed"
          ) %>%
          pull(
            all_of(metric)
          )
      ),
      list(
        Family = "Half_precipitation",
        Effect = "Overall_Half",
        Target = "Half",
        Reference = "Ambient",
        target_values =
          dg %>%
          filter(
            PrecipF == "Half"
          ) %>%
          pull(
            all_of(metric)
          ),
        reference_values =
          dg %>%
          filter(
            PrecipF == "Ambient"
          ) %>%
          pull(
            all_of(metric)
          )
      ),
      list(
        Family = "Double_precipitation",
        Effect = "Overall_Double",
        Target = "Double",
        Reference = "Ambient",
        target_values =
          dg %>%
          filter(
            PrecipF == "Double"
          ) %>%
          pull(
            all_of(metric)
          ),
        reference_values =
          dg %>%
          filter(
            PrecipF == "Ambient"
          ) %>%
          pull(
            all_of(metric)
          )
      ),
      list(
        Family = "Clipping",
        Effect = "Overall_Clipping",
        Target = "Clipped",
        Reference = "Unclipped",
        target_values =
          dg %>%
          filter(
            ClipF == "Clipped"
          ) %>%
          pull(
            all_of(metric)
          ),
        reference_values =
          dg %>%
          filter(
            ClipF == "Unclipped"
          ) %>%
          pull(
            all_of(metric)
          )
      )
    )

    for (dd in definitions) {

      overall_cohen_rows[[idx]] <-
        bind_cols(
          tibble(
            Group = g,
            Metric = metric,
            Family =
              dd$Family,
            Effect =
              dd$Effect,
            Target =
              dd$Target,
            Reference =
              dd$Reference
          ),
          standardized_difference(
            dd$target_values,
            dd$reference_values
          )
        )

      idx <- idx + 1
    }
  }
}


overall_cohen <- bind_rows(
  overall_cohen_rows
) %>%
  mutate(
    Fig2_treatment = dplyr::recode(
      Effect,
      "Overall_Warming" = "W",
      "Overall_Half" = "H",
      "Overall_Double" = "D",
      "Overall_Clipping" = "C"
    ),
    Averaged_over = dplyr::case_when(
      Family == "Warming" ~ "Precipitation x Clipping",
      Family == "Half_precipitation" ~ "Warming x Clipping",
      Family == "Double_precipitation" ~ "Warming x Clipping",
      Family == "Clipping" ~ "Warming x Precipitation",
      TRUE ~ NA_character_
    ),
    Raw_marginal_difference =
      Raw_difference,
    Cohen_d_overall =
      Cohen_d_pooled,
    Hedges_g_overall =
      Hedges_g
  )


write_out(
  overall_cohen,
  "07_pTTPL_OVERALL_Cohen_d_Hedges_g.csv"
)


# ------------------------------------------------------------------------------
# 11C. FINAL FIG. 2a OVERALL EFFECT TABLE
# ------------------------------------------------------------------------------

fig2_overall <- overall_cohen %>%
  left_join(
    primary_overall %>%
      select(
        Group,
        Metric,
        Family,
        Contrast,
        Estimate,
        SE,
        df,
        CI_low,
        CI_high,
        t_ratio,
        P_raw,
        P_Holm_four_overall,
        Significance_raw,
        Significance_Holm
      ),
    by = c(
      "Group",
      "Metric",
      "Family"
    )
  ) %>%
  mutate(
    LMM_minus_raw =
      Estimate -
      Raw_difference,
    Same_direction =
      dplyr::case_when(
        is.na(Estimate) |
          is.na(Raw_difference) ~ NA,
        Estimate == 0 |
          Raw_difference == 0 ~ TRUE,
        sign(Estimate) ==
          sign(Raw_difference) ~ TRUE,
        TRUE ~ FALSE
      ),
    Effect_order =
      factor(
        Effect,
        levels = c(
          "Overall_Warming",
          "Overall_Half",
          "Overall_Double",
          "Overall_Clipping"
        )
      )
  ) %>%
  arrange(
    Metric,
    factor(
      Group,
      levels = GROUP_ORDER
    ),
    Effect_order
  )


write_out(
  fig2_overall,
  "08_pTTPL_FIG2a_OVERALL_W_H_D_C_effectsize_and_LMM.csv"
)


fig2_overall_beta_compact <- fig2_overall %>%
  filter(
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  select(
    Group,
    Metric,
    Family,
    Effect,
    Target,
    Reference,
    N_target,
    N_reference,
    Mean_target,
    Mean_reference,
    Raw_difference,
    Raw_marginal_difference,
    Cohen_d_pooled,
    Cohen_d_overall,
    Hedges_g,
    Hedges_g_overall,
    Fig2_treatment,
    Averaged_over,
    Estimate,
    SE,
    df,
    CI_low,
    CI_high,
    P_raw,
    P_Holm_four_overall,
    Significance_Holm,
    LMM_minus_raw,
    Same_direction
  )


write_out(
  fig2_overall_beta_compact,
  "08_pTTPL_FIG2a_OVERALL_beta_compact.csv"
)

# ------------------------------------------------------------------------------
# 12. BLOCK-MATCHED SENSITIVITY
# ------------------------------------------------------------------------------

get_matched_pair <- function(
  d,
  metric,
  target,
  reference
) {

  a <- d %>%
    filter(
      Treatment == target
    ) %>%
    transmute(
      Block,
      Target_value =
        .data[[metric]],
      WholePlot_target =
        as.character(
          WholePlotID
        ),
      Plot_target =
        as.character(Plot)
    )

  b <- d %>%
    filter(
      Treatment == reference
    ) %>%
    transmute(
      Block,
      Reference_value =
        .data[[metric]],
      WholePlot_reference =
        as.character(
          WholePlotID
        ),
      Plot_reference =
        as.character(Plot)
    )

  inner_join(
    a,
    b,
    by = "Block"
  ) %>%
    mutate(
      Difference =
        Target_value -
        Reference_value,
      SameWholePlot =
        WholePlot_target ==
        WholePlot_reference
    )
}

safe_t_onesample <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (
    length(x) < 2 ||
      !is.finite(sd(x)) ||
      sd(x) == 0
  ) {
    return(list(
      estimate =
        ifelse(
          length(x) > 0,
          mean(x),
          NA_real_
        ),
      se = NA_real_,
      df = length(x) - 1,
      lo = NA_real_,
      hi = NA_real_,
      p = NA_real_
    ))
  }

  tt <- t.test(
    x,
    mu = 0
  )

  list(
    estimate = mean(x),
    se =
      sd(x) /
      sqrt(length(x)),
    df =
      unname(
        tt$parameter
      ),
    lo =
      tt$conf.int[1],
    hi =
      tt$conf.int[2],
    p =
      tt$p.value
  )
}

block_rows <- list()
idx <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    )

  for (metric in MAIN_METRICS) {

    for (
      j in seq_len(
        nrow(contrast_def)
      )
    ) {

      cd <- contrast_def[j, ]

      pp <- get_matched_pair(
        dg,
        metric,
        cd$Target,
        cd$Reference
      )

      tt <- safe_t_onesample(
        pp$Difference
      )

      dz <- if (
        nrow(pp) >= 2 &&
        is.finite(
          sd(pp$Difference)
        ) &&
        sd(pp$Difference) > 0
      ) {
        mean(
          pp$Difference
        ) /
          sd(
            pp$Difference
          )
      } else {
        NA_real_
      }

      block_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        Family =
          cd$Family,
        Contrast_order =
          cd$Contrast_order,
        Contrast =
          cd$Contrast,
        Target =
          cd$Target,
        Reference =
          cd$Reference,
        Context =
          cd$Context,
        N_blocks =
          nrow(pp),
        Mean_block_difference =
          tt$estimate,
        SE_block_difference =
          tt$se,
        df =
          tt$df,
        CI_low =
          tt$lo,
        CI_high =
          tt$hi,
        P_blockmatched =
          tt$p,
        Cohen_dz_blockmatched =
          dz,
        All_same_wholeplot =
          if (
            nrow(pp) > 0
          ) {
            all(
              pp$SameWholePlot
            )
          } else {
            NA
          }
      )

      idx <- idx + 1
    }
  }
}

blockmatched <- bind_rows(
  block_rows
) %>%
  group_by(
    Group,
    Metric,
    Family
  ) %>%
  mutate(
    P_blockmatched_Holm =
      p.adjust(
        P_blockmatched,
        method = "holm"
      )
  ) %>%
  ungroup()

write_out(
  blockmatched,
  "08_pTTPL_blockmatched_effect_specific_sensitivity.csv"
)

# ------------------------------------------------------------------------------
# 13. BETA_p SLOPE-SE PROPAGATION
# ------------------------------------------------------------------------------
# The uncertainty of each subplot-level beta_p estimate is propagated through
# V = SE_Beta^2 in rma.mv(). The fixed-effect structure mirrors the categorical
# primary LMM so Half and Double precipitation are estimated independently.
# This remains a SENSITIVITY ANALYSIS and does not replace the primary LMM.
# ------------------------------------------------------------------------------

fit_beta_se_model <- function(
  d,
  beta_col,
  se_col
) {

  dd <- d %>%
    filter(
      is.finite(
        .data[[beta_col]]
      ),
      is.finite(
        .data[[se_col]]
      ),
      .data[[se_col]] > 0
    ) %>%
    droplevels()

  if (nrow(dd) < 20) {
    return(list(
      fit = NULL,
      data = dd,
      error =
        "Too few observations."
    ))
  }

  fit <- tryCatch(
    metafor::rma.mv(
      yi =
        dd[[beta_col]],
      V =
        dd[[se_col]]^2,
      mods =
        ~ WarmF *
        PrecipF *
        ClipF,
      random = list(
        ~ 1 | Block,
        ~ 1 | BlockWholePlotID,
        ~ 1 | Plot
      ),
      data = dd,
      method = "REML",
      test = "t",
      dfs = "contain"
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      fit = NULL,
      data = dd,
      error =
        conditionMessage(fit)
    ))
  }

  list(
    fit = fit,
    data = dd,
    error = NA_character_
  )
}


# One-row factor representation of each manuscript treatment label.
# Levels are fixed explicitly so model.matrix() reproduces the primary
# WarmF * PrecipF * ClipF parameterization.
treatment_factor_row <- function(
  label
) {

  grid <- tibble(
    Treatment =
      TREATMENT_LEVELS,
    WarmF = factor(
      c(
        "Unwarmed", "Warmed", "Unwarmed", "Unwarmed", "Unwarmed",
        "Warmed", "Warmed", "Warmed", "Unwarmed", "Unwarmed",
        "Warmed", "Warmed"
      ),
      levels = c(
        "Unwarmed",
        "Warmed"
      )
    ),
    PrecipF = factor(
      c(
        "Ambient", "Ambient", "Half", "Double", "Ambient",
        "Half", "Double", "Ambient", "Half", "Double",
        "Half", "Double"
      ),
      levels = c(
        "Ambient",
        "Half",
        "Double"
      )
    ),
    ClipF = factor(
      c(
        "Unclipped", "Unclipped", "Unclipped", "Unclipped", "Clipped",
        "Unclipped", "Unclipped", "Clipped", "Clipped", "Clipped",
        "Clipped", "Clipped"
      ),
      levels = c(
        "Unclipped",
        "Clipped"
      )
    )
  )

  out <- grid %>%
    filter(
      Treatment ==
        label
    )

  if (nrow(out) != 1) {
    stop(
      "Unknown treatment label: ",
      label
    )
  }

  out
}


meta_pair_weight <- function(
  fit,
  target,
  reference
) {

  tt <- treatment_factor_row(
    target
  )

  rr <- treatment_factor_row(
    reference
  )

  mm <- model.matrix(
    ~ WarmF *
      PrecipF *
      ClipF,
    data = bind_rows(
      tt,
      rr
    )
  )

  drow <-
    mm[1, , drop = TRUE] -
    mm[2, , drop = TRUE]

  bn <- names(
    coef(fit)
  )

  L <- rep(
    0,
    length(bn)
  )

  names(L) <- bn

  for (nm in bn) {

    source_name <- if (
      nm == "intrcpt"
    ) {
      "(Intercept)"
    } else {
      nm
    }

    if (
      source_name %in%
        names(drow)
    ) {

      L[nm] <-
        drow[source_name]

    } else {

      stop(
        "Could not map meta-regression coefficient '",
        nm,
        "' to the categorical WarmF*PrecipF*ClipF design matrix."
      )
    }
  }

  L
}


meta_model_contrast <- function(
  fit,
  L,
  label
) {

  Lmat <- matrix(
    as.numeric(L),
    nrow = 1,
    dimnames = list(
      label,
      names(L)
    )
  )

  aa <- tryCatch(
    stats::anova(
      fit,
      X = Lmat
    ),
    error = function(e) e
  )

  bb <- tryCatch(
    as.numeric(
      stats::coef(fit)
    ),
    error = function(e)
      NULL
  )

  vv <- tryCatch(
    as.matrix(
      stats::vcov(fit)
    ),
    error = function(e)
      NULL
  )

  est <- NA_real_
  se <- NA_real_

  if (
    !is.null(bb) &&
      length(bb) ==
      length(L)
  ) {
    est <- sum(
      as.numeric(L) *
        bb
    )
  }

  if (
    !is.null(vv) &&
      nrow(vv) ==
      length(L) &&
      ncol(vv) ==
      length(L)
  ) {

    se2 <- as.numeric(
      t(as.numeric(L)) %*%
        vv %*%
        as.numeric(L)
    )

    if (
      is.finite(se2) &&
      se2 >= 0
    ) {
      se <- sqrt(se2)
    }
  }

  p <- NA_real_
  ddf <- NA_real_
  test_error <- NA_character_

  if (inherits(aa, "error")) {

    test_error <-
      conditionMessage(aa)

  } else {

    if (!is.null(aa$pval)) {
      p <- as.numeric(
        aa$pval[1]
      )
    } else if (
      !is.null(aa$QMp)
    ) {
      p <- as.numeric(
        aa$QMp[1]
      )
    }

    if (!is.null(aa$ddf)) {
      ddf <- as.numeric(
        aa$ddf[1]
      )
    } else if (
      !is.null(aa$dfs)
    ) {
      ddf <- as.numeric(
        aa$dfs[1]
      )
    } else if (
      !is.null(aa$QMdf) &&
      length(aa$QMdf) >= 2
    ) {
      ddf <- as.numeric(
        aa$QMdf[2]
      )
    }
  }

  if (
    !is.finite(p) &&
    is.finite(est) &&
    is.finite(se) &&
    se > 0
  ) {

    stat <- est / se

    if (
      is.finite(ddf) &&
      ddf > 0
    ) {
      p <- 2 *
        stats::pt(
          abs(stat),
          df = ddf,
          lower.tail = FALSE
        )
    } else {
      p <- 2 *
        stats::pnorm(
          abs(stat),
          lower.tail = FALSE
        )
    }
  }

  crit <- if (
    is.finite(ddf) &&
    ddf > 0
  ) {
    stats::qt(
      0.975,
      df = ddf
    )
  } else {
    stats::qnorm(
      0.975
    )
  }

  lo <- if (
    is.finite(est) &&
    is.finite(se)
  ) {
    est - crit * se
  } else {
    NA_real_
  }

  hi <- if (
    is.finite(est) &&
    is.finite(se)
  ) {
    est + crit * se
  } else {
    NA_real_
  }

  tibble(
    Contrast = label,
    Estimate_SEprop = est,
    SE_SEprop = se,
    CI_low_SEprop = lo,
    CI_high_SEprop = hi,
    df_SEprop = ddf,
    P_SEprop = p,
    Test_error = test_error
  )
}

meta_diag_rows <- list()
meta_contrast_rows <- list()
meta_overall_rows <- list()
meta_coef_rows <- list()

idx_diag <- 1
idx_con <- 1
idx_over <- 1
idx_coef <- 1

for (g in GROUP_ORDER) {

  dg <- dat %>%
    filter(
      Group == g
    )

  for (
    beta_metric in names(
      BETA_SE_PAIRS
    )
  ) {

    se_metric <- unname(
      BETA_SE_PAIRS[
        beta_metric
      ]
    )

    ans <- fit_beta_se_model(
      dg,
      beta_metric,
      se_metric
    )

    if (is.null(ans$fit)) {

      meta_diag_rows[[idx_diag]] <-
        tibble(
          Group = g,
          Metric = beta_metric,
          SE_metric = se_metric,
          N = nrow(ans$data),
          N_blocks =
            n_distinct(
              ans$data$Block
            ),
          N_wholeplots =
            n_distinct(
              ans$data$
                BlockWholePlotID
            ),
          Median_input_SE =
            if (
              nrow(ans$data) > 0
            ) {
              median(
                ans$data[[
                  se_metric
                ]],
                na.rm = TRUE
              )
            } else {
              NA_real_
            },
          Mean_input_SE =
            if (
              nrow(ans$data) > 0
            ) {
              mean(
                ans$data[[
                  se_metric
                ]],
                na.rm = TRUE
              )
            } else {
              NA_real_
            },
          Empirical_SD_Beta =
            if (
              nrow(ans$data) > 1
            ) {
              sd(
                ans$data[[
                  beta_metric
                ]],
                na.rm = TRUE
              )
            } else {
              NA_real_
            },
          Block_variance = NA_real_,
          WholePlot_within_Block_variance =
            NA_real_,
          Plot_heterogeneity_variance =
            NA_real_,
          QE = NA_real_,
          QE_P = NA_real_,
          Inference_method =
            "rma.mv; V=SE_Beta^2; categorical treatment factors; test=t; dfs=contain",
          Error = ans$error
        )

      idx_diag <- idx_diag + 1
      next
    }

    fit <- ans$fit
    dd <- ans$data
    sig2 <- fit$sigma2

    meta_diag_rows[[idx_diag]] <-
      tibble(
        Group = g,
        Metric = beta_metric,
        SE_metric = se_metric,
        N = nrow(dd),
        N_blocks =
          n_distinct(dd$Block),
        N_wholeplots =
          n_distinct(
            dd$BlockWholePlotID
          ),
        Median_input_SE =
          median(
            dd[[se_metric]],
            na.rm = TRUE
          ),
        Mean_input_SE =
          mean(
            dd[[se_metric]],
            na.rm = TRUE
          ),
        Empirical_SD_Beta =
          sd(
            dd[[beta_metric]],
            na.rm = TRUE
          ),
        Block_variance =
          ifelse(
            length(sig2) >= 1,
            sig2[1],
            NA
          ),
        WholePlot_within_Block_variance =
          ifelse(
            length(sig2) >= 2,
            sig2[2],
            NA
          ),
        Plot_heterogeneity_variance =
          ifelse(
            length(sig2) >= 3,
            sig2[3],
            NA
          ),
        QE = fit$QE,
        QE_P = fit$QEp,
        Inference_method =
          "rma.mv; V=SE_Beta^2; categorical treatment factors; test=t; dfs=contain",
        Error = NA_character_
      )

    idx_diag <- idx_diag + 1

    cf <- tibble(
      Group = g,
      Metric = beta_metric,
      Coefficient =
        names(coef(fit)),
      Estimate =
        as.numeric(
          coef(fit)
        ),
      SE =
        as.numeric(
          fit$se
        ),
      t_or_z =
        as.numeric(
          fit$zval
        ),
      P =
        as.numeric(
          fit$pval
        ),
      CI_low =
        as.numeric(
          fit$ci.lb
        ),
      CI_high =
        as.numeric(
          fit$ci.ub
        )
    )

    meta_coef_rows[[idx_coef]] <-
      cf
    idx_coef <- idx_coef + 1

    family_weights <- list()

    for (
      fam in unique(
        contrast_def$Family
      )
    ) {

      defs <- contrast_def %>%
        filter(
          Family == fam
        ) %>%
        arrange(
          Contrast_order
        )

      fw <- list()

      for (
        j in seq_len(
          nrow(defs)
        )
      ) {

        L <- meta_pair_weight(
          fit,
          defs$Target[j],
          defs$Reference[j]
        )

        fw[[j]] <- L

        rr <- meta_model_contrast(
          fit,
          L,
          defs$Contrast[j]
        ) %>%
          mutate(
            Group = g,
            Metric =
              beta_metric,
            Family = fam,
            Contrast_order =
              defs$Contrast_order[j],
            Target =
              defs$Target[j],
            Reference =
              defs$Reference[j],
            Context =
              defs$Context[j],
            .before = 1
          )

        meta_contrast_rows[[idx_con]] <-
          rr
        idx_con <- idx_con + 1
      }

      family_weights[[fam]] <- fw

      Lavg <- make_average_weight(
        fw
      )

      oo <- meta_model_contrast(
        fit,
        Lavg,
        paste0(
          "Overall_",
          fam
        )
      ) %>%
        mutate(
          Group = g,
          Metric =
            beta_metric,
          Family = fam,
          .before = 1
        )

      meta_overall_rows[[idx_over]] <-
        oo
      idx_over <- idx_over + 1
    }
  }
}

meta_diag <- bind_rows(
  meta_diag_rows
)

meta_coefficients <- bind_rows(
  meta_coef_rows
)

meta_contrasts <- bind_rows(
  meta_contrast_rows
) %>%
  group_by(
    Group,
    Metric,
    Family
  ) %>%
  mutate(
    P_SEprop_Holm_within_family =
      p.adjust(
        P_SEprop,
        method = "holm"
      ),
    Significance_SEprop_Holm =
      sig_label(
        P_SEprop_Holm_within_family
      )
  ) %>%
  ungroup() %>%
  group_by(
    Group,
    Metric
  ) %>%
  mutate(
    P_SEprop_Holm_all20 =
      p.adjust(
        P_SEprop,
        method = "holm"
      )
  ) %>%
  ungroup()

meta_overall <- bind_rows(
  meta_overall_rows
) %>%
  group_by(
    Group,
    Metric
  ) %>%
  mutate(
    P_SEprop_Holm_four_overall =
      p.adjust(
        P_SEprop,
        method = "holm"
      ),
    Significance_SEprop_raw =
      sig_label(P_SEprop)
  ) %>%
  ungroup()

write_out(
  meta_diag,
  "10_pTTPL_Beta_SE_propagation_model_diagnostics.csv"
)

write_out(
  meta_coefficients,
  "10_pTTPL_Beta_SE_propagation_coefficients.csv"
)

write_out(
  meta_contrasts,
  "10_pTTPL_Beta_SE_propagation_effect_specific_contrasts.csv"
)

write_out(
  meta_overall,
  "10_pTTPL_Beta_SE_propagation_overall_directional_effects.csv"
)

writeLines(
  c(
    "CR2 was intentionally not used.",
    "beta_p estimation uncertainty is propagated through V = SE_Beta^2 in rma.mv().",
    "The meta-regression uses the same Warming/Precipitation/Clipping fixed-effect structure",
    "and random effects for Block, whole plot within Block, and plot-level heterogeneity.",
    "With only four top-level Blocks, a block-clustered CR2 small-sample analysis is not appropriate."
  ),
  file.path(
    OUTDIR,
    "10_pTTPL_Beta_SE_CR2_NOT_USED.txt"
  )
)

# ------------------------------------------------------------------------------
# 14. UNIFIED TABLES
# ------------------------------------------------------------------------------

unified_effect_specific <- primary_contrasts %>%
  left_join(
    cohen_out %>%
      select(
        Group,
        Metric,
        Family,
        Contrast,
        Cohen_d_pooled,
        Hedges_g,
        Mean_target,
        Mean_reference
      ),
    by = c(
      "Group",
      "Metric",
      "Family",
      "Contrast"
    )
  ) %>%
  left_join(
    blockmatched %>%
      select(
        Group,
        Metric,
        Family,
        Contrast,
        Mean_block_difference,
        SE_block_difference,
        P_blockmatched,
        P_blockmatched_Holm,
        Cohen_dz_blockmatched,
        All_same_wholeplot
      ),
    by = c(
      "Group",
      "Metric",
      "Family",
      "Contrast"
    )
  )

if (nrow(meta_contrasts) > 0) {

  unified_effect_specific <-
    unified_effect_specific %>%
    left_join(
      meta_contrasts %>%
        select(
          Group,
          Metric,
          Family,
          Contrast,
          Estimate_SEprop,
          SE_SEprop,
          CI_low_SEprop,
          CI_high_SEprop,
          df_SEprop,
          P_SEprop,
          P_SEprop_Holm_within_family
        ),
      by = c(
        "Group",
        "Metric",
        "Family",
        "Contrast"
      )
    )
}

write_out(
  unified_effect_specific,
  "11_pTTPL_UNIFIED_effect_specific_all_evidence.csv"
)

unified_overall <- primary_overall %>%
  left_join(
    overall_cohen %>%
      select(
        Group,
        Metric,
        Family,
        Effect,
        Target,
        Reference,
        N_target,
        N_reference,
        Mean_target,
        Mean_reference,
        Raw_difference,
        Cohen_d_pooled,
        Hedges_g
      ),
    by = c(
      "Group",
      "Metric",
      "Family"
    )
  )

if (nrow(meta_overall) > 0) {

  unified_overall <-
    unified_overall %>%
    left_join(
      meta_overall %>%
        select(
          Group,
          Metric,
          Family,
          Contrast,
          Estimate_SEprop,
          SE_SEprop,
          CI_low_SEprop,
          CI_high_SEprop,
          df_SEprop,
          P_SEprop,
          P_SEprop_Holm_four_overall
        ),
      by = c(
        "Group",
        "Metric",
        "Family",
        "Contrast"
      )
    )
}

write_out(
  unified_overall,
  "11_pTTPL_UNIFIED_overall_directional_effects.csv"
)

unified_vsctrl <- vsctrl %>%
  left_join(
    cohen_vsctrl %>%
      select(
        Group,
        Metric,
        Treatment,
        Contrast,
        Cohen_d_pooled,
        Hedges_g,
        Mean_target,
        Mean_reference
      ),
    by = c(
      "Group",
      "Metric",
      "Treatment",
      "Contrast"
    )
  )

write_out(
  unified_vsctrl,
  "11_pTTPL_UNIFIED_all_treatments_vs_Ctrl.csv"
)

# ------------------------------------------------------------------------------
# 15. beta_p / REFERENCE-VARIANCE CO-VARIATION UNDER WARMING
# ------------------------------------------------------------------------------

run_covary_analysis <- function(
  d,
  beta_col,
  logvar_col,
  label
) {

  required <- c(
    beta_col,
    logvar_col
  )

  if (!all(
    required %in% names(d)
  )) {
    return(list(
      data = tibble(),
      results = tibble()
    ))
  }

  unwarmed <- d %>%
    filter(
      Warming == 0,
      is.finite(
        .data[[beta_col]]
      ),
      is.finite(
        .data[[logvar_col]]
      )
    ) %>%
    transmute(
      Group,
      Block,
      PrecipF,
      ClipF,
      Beta_unwarmed =
        .data[[beta_col]],
      LogVar_ref_unwarmed =
        .data[[logvar_col]],
      WholePlot_unwarmed =
        as.character(
          WholePlotID
        ),
      Plot_unwarmed =
        as.character(Plot)
    )

  warmed <- d %>%
    filter(
      Warming == 1,
      is.finite(
        .data[[beta_col]]
      ),
      is.finite(
        .data[[logvar_col]]
      )
    ) %>%
    transmute(
      Group,
      Block,
      PrecipF,
      ClipF,
      Beta_warmed =
        .data[[beta_col]],
      LogVar_ref_warmed =
        .data[[logvar_col]],
      WholePlot_warmed =
        as.character(
          WholePlotID
        ),
      Plot_warmed =
        as.character(Plot)
    )

  delta <- inner_join(
    warmed,
    unwarmed,
    by = c(
      "Group",
      "Block",
      "PrecipF",
      "ClipF"
    )
  ) %>%
    mutate(
      Delta_Beta =
        Beta_warmed -
        Beta_unwarmed,
      Delta_LogVar_ref =
        LogVar_ref_warmed -
        LogVar_ref_unwarmed,
      PairID = interaction(
        Block,
        PrecipF,
        drop = TRUE,
        sep = ":P"
      ),
      Context = paste0(
        as.character(
          PrecipF
        ),
        " | ",
        as.character(
          ClipF
        )
      ),
      Variant = label
    ) %>%
    arrange(
      Group,
      Block,
      PrecipF,
      ClipF
    )

  result_rows <- list()
  ii <- 1

  for (g in GROUP_ORDER) {

    dd <- delta %>%
      filter(
        Group == g,
        is.finite(
          Delta_Beta
        ),
        is.finite(
          Delta_LogVar_ref
        )
      ) %>%
      mutate(
        Block =
          droplevels(
            factor(Block)
          ),
        PairID =
          droplevels(
            factor(PairID)
          )
      )

    if (
      nrow(dd) < 12 ||
      !is.finite(
        sd(dd$Delta_Beta)
      ) ||
      sd(dd$Delta_Beta) == 0
    ) {

      result_rows[[ii]] <- tibble(
        Variant = label,
        Group = g,
        N = nrow(dd),
        N_blocks =
          n_distinct(dd$Block),
        N_pairIDs =
          n_distinct(dd$PairID),
        Model =
          "Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID)",
        Null_hypothesis =
          "H0: coefficient(Delta_Beta) = 0",
        Estimate_gamma = NA_real_,
        SE = NA_real_,
        df = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        t_ratio = NA_real_,
        P_raw = NA_real_,
        PairID_variance =
          NA_real_,
        Residual_variance =
          NA_real_,
        Singular = NA,
        OLS_Estimate_gamma =
          NA_real_,
        OLS_SE = NA_real_,
        OLS_df = NA_real_,
        OLS_t = NA_real_,
        OLS_P = NA_real_,
        Error =
          "Too few observations or zero SD in Delta_Beta."
      )

      ii <- ii + 1
      next
    }

    mixed <- tryCatch(
      lmerTest::lmer(
        Delta_LogVar_ref ~
          Delta_Beta +
          PrecipF *
          ClipF +
          Block +
          (1 | PairID),
        data = dd,
        REML = TRUE
      ),
      error = function(e) e
    )

    ols <- tryCatch(
      lm(
        Delta_LogVar_ref ~
          Delta_Beta +
          PrecipF *
          ClipF +
          Block,
        data = dd
      ),
      error = function(e) e
    )

    if (inherits(mixed, "error")) {

      est <- se <- ddf <- tv <- pv <-
        lo <- hi <- pair_var <-
        resid_var <- NA_real_

      singular <- NA
      mixed_error <-
        conditionMessage(mixed)

    } else {

      sm <- as.data.frame(
        summary(mixed)$coefficients
      )

      rr <- sm[
        "Delta_Beta",
        ,
        drop = FALSE
      ]

      est <-
        rr[1, "Estimate"]
      se <-
        rr[1, "Std. Error"]
      ddf <-
        rr[1, "df"]
      tv <-
        rr[1, "t value"]
      pv <-
        rr[1, "Pr(>|t|)"]

      crit <- qt(
        0.975,
        df = ddf
      )

      lo <- est - crit * se
      hi <- est + crit * se

      vc <- as.data.frame(
        VarCorr(mixed)
      )

      pair_var <- vc %>%
        filter(
          grp == "PairID"
        ) %>%
        pull(vcov)

      if (length(pair_var) == 0) {
        pair_var <- NA_real_
      }

      resid_var <- vc %>%
        filter(
          grp == "Residual"
        ) %>%
        pull(vcov)

      if (length(resid_var) == 0) {
        resid_var <- NA_real_
      }

      singular <-
        lme4::isSingular(
          mixed,
          tol = 1e-6
        )

      mixed_error <-
        NA_character_
    }

    if (inherits(ols, "error")) {

      oe <- os <- odf <- ot <- op <-
        NA_real_

    } else {

      so <- summary(ols)
      co <- as.data.frame(
        so$coefficients
      )

      ro <- co[
        "Delta_Beta",
        ,
        drop = FALSE
      ]

      oe <-
        ro[1, "Estimate"]
      os <-
        ro[1, "Std. Error"]
      ot <-
        ro[1, "t value"]
      op <-
        ro[1, "Pr(>|t|)"]
      odf <-
        df.residual(ols)
    }

    result_rows[[ii]] <- tibble(
      Variant = label,
      Group = g,
      N = nrow(dd),
      N_blocks =
        n_distinct(dd$Block),
      N_pairIDs =
        n_distinct(dd$PairID),
      Model =
        "Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID)",
      Null_hypothesis =
        "H0: coefficient(Delta_Beta) = 0",
      Estimate_gamma = est,
      SE = se,
      df = ddf,
      CI_low = lo,
      CI_high = hi,
      t_ratio = tv,
      P_raw = pv,
      PairID_variance =
        pair_var[1],
      Residual_variance =
        resid_var[1],
      Singular = singular,
      OLS_Estimate_gamma = oe,
      OLS_SE = os,
      OLS_df = odf,
      OLS_t = ot,
      OLS_P = op,
      Error = mixed_error
    )

    ii <- ii + 1
  }

  results <- bind_rows(
    result_rows
  ) %>%
    mutate(
      P_Holm_across_groups =
        p.adjust(
          P_raw,
          method = "holm"
        ),
      Significance_Holm =
        sig_label(
          P_Holm_across_groups
        )
    )

  list(
    data = delta,
    results = results
  )
}

cov_raw <- run_covary_analysis(
  dat,
  "Beta_raw",
  "LogVar_at_ref_raw",
  "Raw"
)

cov_det <- run_covary_analysis(
  dat,
  "Beta_detrended",
  "LogVar_at_ref_detrended",
  "Detrended"
)

covary_data <- bind_rows(
  cov_raw$data,
  cov_det$data
)

covary_results <- bind_rows(
  cov_raw$results,
  cov_det$results
)

write_out(
  covary_data,
  "12_pTTPL_Beta_LogVar_covary_warming_change_data.csv"
)

write_out(
  covary_results,
  "12_pTTPL_Beta_LogVar_covary_adjusted_tests.csv"
)

# ------------------------------------------------------------------------------
# 16. OPTIONAL beta_p — WHOLE-COMMUNITY TEMPORAL VARIABILITY ASSOCIATION
# ------------------------------------------------------------------------------

ctv_assoc <- tibble()

if (
  RUN_CTV_ASSOCIATION &&
  file.exists(CTV_SOURCE_FILE)
) {

  ctv <- read.csv(
    CTV_SOURCE_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_ctv <- c(
    "Group",
    "Plot",
    "CTV_Hellinger_RMS"
  )

  if (
    all(
      required_ctv %in%
        names(ctv)
    )
  ) {

    ctv_small <- ctv %>%
      transmute(
        Group =
          as.character(Group),
        Plot =
          as.character(Plot),
        CTV_Hellinger_RMS =
          as.numeric(
            CTV_Hellinger_RMS
          )
      )

    dd <- dat %>%
      mutate(
        Group_chr =
          as.character(Group),
        Plot_chr =
          as.character(Plot)
      ) %>%
      left_join(
        ctv_small,
        by = c(
          "Group_chr" =
            "Group",
          "Plot_chr" =
            "Plot"
        )
      )

    assoc_rows <- list()
    ii <- 1

    for (g in GROUP_ORDER) {

      dg <- dd %>%
        filter(
          Group_chr == g,
          is.finite(
            Beta_raw
          ),
          is.finite(
            CTV_Hellinger_RMS
          )
        ) %>%
        mutate(
          Beta_p_raw_z =
            as.numeric(
              scale(Beta_raw)
            )
        ) %>%
        droplevels()

      m <- tryCatch(
        lmerTest::lmer(
          CTV_Hellinger_RMS ~
            Beta_p_raw_z +
            WarmF *
            PrecipF *
            ClipF +
            (1 | Block/WholePlotID),
          data = dg,
          REML = TRUE
        ),
        error = function(e) e
      )

      if (!inherits(m, "error")) {

        co <- as.data.frame(
          summary(m)$coefficients
        )

        rr <- co[
          "Beta_p_raw_z",
          ,
          drop = FALSE
        ]

        crit <- qt(
          0.975,
          df =
            rr[1, "df"]
        )

        assoc_rows[[ii]] <-
          tibble(
            Group = g,
            N = nrow(dg),
            Estimate_CTV_per_1SD_Beta_p =
              rr[1, "Estimate"],
            SE =
              rr[1, "Std. Error"],
            df =
              rr[1, "df"],
            CI_low =
              rr[1, "Estimate"] -
              crit *
              rr[1, "Std. Error"],
            CI_high =
              rr[1, "Estimate"] +
              crit *
              rr[1, "Std. Error"],
            t =
              rr[1, "t value"],
            P =
              rr[1, "Pr(>|t|)"],
            Singular =
              lme4::isSingular(
                m,
                tol = 1e-6
              )
          )

        ii <- ii + 1
      }
    }

    ctv_assoc <- bind_rows(
      assoc_rows
    ) %>%
      mutate(
        P_Holm_across_groups =
          p.adjust(
            P,
            method = "holm"
          )
      )

    write_out(
      ctv_assoc,
      "13_pTTPL_beta_CTV_association_adjusted_descriptive.csv"
    )

  } else {

    warning(
      "CTV source exists but is missing required columns; ",
      "beta_p–CTV association was skipped."
    )
  }

} else if (RUN_CTV_ASSOCIATION) {

  warning(
    "CTV source file not found; beta_p–CTV association was skipped: ",
    CTV_SOURCE_FILE
  )
}

# ------------------------------------------------------------------------------
# 17. COMPACT TABLES FOR MANUSCRIPT / RESPONSE
# ------------------------------------------------------------------------------

overall_warming_compact <-
  unified_overall %>%
  filter(
    Family == "Warming",
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  arrange(
    Metric,
    factor(
      Group,
      levels =
        GROUP_ORDER
    )
  )

write_out(
  overall_warming_compact,
  "14_pTTPL_REVIEWER_overall_warming_compact.csv"
)

reference_variance_warming <-
  primary_overall %>%
  filter(
    Family == "Warming",
    Metric %in% c(
      "LogVar_at_ref_raw",
      "LogVar_at_ref_detrended"
    )
  ) %>%
  arrange(
    Metric,
    factor(
      Group,
      levels =
        GROUP_ORDER
    )
  )

write_out(
  reference_variance_warming,
  "14_pTTPL_REVIEWER_reference_variance_warming.csv"
)

beta_context_compact <-
  unified_effect_specific %>%
  filter(
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  arrange(
    Metric,
    Family,
    factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Contrast_order
  )

write_out(
  beta_context_compact,
  "14_pTTPL_REVIEWER_beta_effect_specific_all_evidence.csv"
)

# Fig. 2a uses OVERALL factor effects, not the individual W/H/D/C-vs-Ctrl
# treatment contrasts. This table aligns pTTPL directly with the revised tTTPL
# Fig. 2a analysis.
fig2a_overall_inference <- fig2_overall %>%
  filter(
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  arrange(
    Metric,
    factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Effect_order
  )

write_out(
  fig2a_overall_inference,
  "14_pTTPL_REVIEWER_Fig2a_OVERALL_effects.csv"
)


# Table 1 remains individual treatment-combination vs Ctrl. Do NOT use the
# Fig. 2a overall P values for these rows.
table1_inference <- unified_vsctrl %>%
  filter(
    Metric == "Beta_raw"
  ) %>%
  arrange(
    factor(
      Treatment,
      levels =
        setdiff(
          TREATMENT_LEVELS,
          "Ctrl"
        )
    ),
    factor(
      Group,
      levels =
        GROUP_ORDER
    )
  )

write_out(
  table1_inference,
  "14_pTTPL_REVIEWER_Table1_Fig2_vsCtrl_inference.csv"
)

# ------------------------------------------------------------------------------
# 18. FIGURES
# ------------------------------------------------------------------------------

# 18A. 12-treatment pTTPL beta overview.
dplot <- descriptive %>%
  filter(
    Metric == "Beta_raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Treatment = factor(
      Treatment,
      levels =
        TREATMENT_LEVELS
    )
  )

if (nrow(dplot) > 0) {

  p <- ggplot(
    dplot,
    aes(
      x = Treatment,
      y = Mean,
      group = 1
    )
  ) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = Mean - SE,
        ymax = Mean + SE
      ),
      width = 0.2
    ) +
    facet_wrap(
      ~ Group,
      scales = "free_y"
    ) +
    theme_bw() +
    labs(
      x = NULL,
      y = expression(beta[p]),
      title =
        "pTTPL exponent across the 12 experimental treatments"
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_QC_12treatments_Beta_raw",
    11,
    5.5
  )
}

# 18B. Fig. 2a overall factor effects for pTTPL beta_p.
#
# Bar magnitude = descriptive OVERALL Cohen's d.
# Asterisks = Holm-adjusted P from the corresponding categorical LMM overall
#             contrast (four overall factors corrected within Group x Metric).
fig2a_plot <- fig2_overall %>%
  filter(
    Metric == "Beta_raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Effect_label = dplyr::recode(
      Effect,
      "Overall_Warming" = "Warming",
      "Overall_Half" = "Half precipitation",
      "Overall_Double" = "Double precipitation",
      "Overall_Clipping" = "Clipping"
    ),
    Effect_label = factor(
      Effect_label,
      levels = c(
        "Warming",
        "Half precipitation",
        "Double precipitation",
        "Clipping"
      )
    ),
    Plot_sig =
      ifelse(
        Significance_Holm %in%
          c(
            "*",
            "**",
            "***"
          ),
        Significance_Holm,
        ""
      )
  )

if (nrow(fig2a_plot) > 0) {

  p <- ggplot(
    fig2a_plot,
    aes(
      x = Effect_label,
      y = Cohen_d_overall
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +
    geom_col(
      width = 0.70
    ) +
    geom_text(
      aes(
        label = Plot_sig
      ),
      nudge_y =
        0.08,
      size = 4
    ) +
    facet_wrap(
      ~ Group,
      nrow = 1
    ) +
    theme_bw() +
    labs(
      x = NULL,
      y = "Overall Cohen's d",
      title =
        "Overall factor effects on pTTPL beta"
    ) +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_pdf_png(
    p,
    "Fig2a_pTTPL_OVERALL_Cohen_d_with_categorical_LMM_Holm",
    12,
    4.8
  )
}


# QC forest plot: categorical-LMM overall Delta beta_p with 95% CI.
fig2a_lmm_plot <- fig2_overall %>%
  filter(
    Metric == "Beta_raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Effect_label = dplyr::recode(
      Effect,
      "Overall_Warming" = "Warming",
      "Overall_Half" = "Half precipitation",
      "Overall_Double" = "Double precipitation",
      "Overall_Clipping" = "Clipping"
    ),
    Effect_label = factor(
      Effect_label,
      levels = rev(
        c(
          "Warming",
          "Half precipitation",
          "Double precipitation",
          "Clipping"
        )
      )
    )
  )

if (nrow(fig2a_lmm_plot) > 0) {

  p <- ggplot(
    fig2a_lmm_plot,
    aes(
      x = Estimate,
      y = Effect_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.18
    ) +
    geom_point() +
    facet_wrap(
      ~ Group,
      nrow = 1,
      scales = "free_x"
    ) +
    theme_bw() +
    labs(
      x =
        expression(
          Delta * beta[p] ~
            "(overall categorical-LMM contrast)"
        ),
      y = NULL
    )

  save_pdf_png(
    p,
    "Fig2a_pTTPL_OVERALL_categorical_LMM_contrasts",
    12,
    4.8
  )
}


# 18B. Overall warming: raw vs detrended beta_p.
pw <- primary_overall %>%
  filter(
    Family == "Warming",
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Metric_label =
      ifelse(
        Metric == "Beta_raw",
        "Original",
        "Detrended"
      )
  )

if (nrow(pw) > 0) {

  p <- ggplot(
    pw,
    aes(
      x = Estimate,
      y = Group,
      shape = Metric_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.2,
      position =
        position_dodge(
          width = 0.5
        )
    ) +
    geom_point(
      position =
        position_dodge(
          width = 0.5
        )
    ) +
    theme_bw() +
    labs(
      x =
        expression(
          Delta * beta[p] ~
            "(warming effect)"
        ),
      y = NULL,
      shape = NULL
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_overall_warming_raw_detrended",
    7,
    4.5
  )
}

# 18C. Reference-level variance warming effect.
pr <- primary_overall %>%
  filter(
    Family == "Warming",
    Metric ==
      "LogVar_at_ref_raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    )
  )

if (nrow(pr) > 0) {

  p <- ggplot(
    pr,
    aes(
      x = Estimate,
      y = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.2
    ) +
    geom_point() +
    theme_bw() +
    labs(
      x =
        "Warming effect on log10 variance at the pTTPL reference level",
      y = NULL
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_reference_variance_overall_warming",
    7,
    4.5
  )
}

# 18D. Co-variation gamma forest plot.
pc <- covary_results %>%
  filter(
    Variant == "Raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    )
  )

if (nrow(pc) > 0) {

  p <- ggplot(
    pc,
    aes(
      x = Estimate_gamma,
      y = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.2
    ) +
    geom_point() +
    theme_bw() +
    labs(
      x =
        expression(
          gamma ~
            "(association between " *
            Delta * beta[p] *
            " and " *
            Delta * logVar[ref] *
            ")"
        ),
      y = NULL
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_Beta_LogVar_covary_gamma",
    7.5,
    4.5
  )
}

# 18E. SE-propagated overall warming.
pse <- meta_overall %>%
  filter(
    Family == "Warming"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels =
        GROUP_ORDER
    ),
    Metric_label =
      ifelse(
        Metric == "Beta_raw",
        "Original",
        "Detrended"
      )
  )

if (nrow(pse) > 0) {

  p <- ggplot(
    pse,
    aes(
      x = Estimate_SEprop,
      y = Group,
      shape = Metric_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low_SEprop,
        xmax = CI_high_SEprop
      ),
      height = 0.2,
      position =
        position_dodge(
          width = 0.5
        )
    ) +
    geom_point(
      position =
        position_dodge(
          width = 0.5
        )
    ) +
    theme_bw() +
    labs(
      x =
        expression(
          Delta * beta[p] ~
            "(SE-propagated warming effect)"
        ),
      y = NULL,
      shape = NULL
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_SEprop_overall_warming",
    7.5,
    4.5
  )
}

# 18F. Optional beta_p–CTV association.
if (nrow(ctv_assoc) > 0) {

  pcv <- ctv_assoc %>%
    mutate(
      Group = factor(
        Group,
        levels =
          GROUP_ORDER
      )
    )

  p <- ggplot(
    pcv,
    aes(
      x =
        Estimate_CTV_per_1SD_Beta_p,
      y = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.2
    ) +
    geom_point() +
    theme_bw() +
    labs(
      x =
        expression(
          Delta * CTV ~
            "per 1 SD increase in " *
            beta[p]
        ),
      y = NULL
    )

  save_pdf_png(
    p,
    "Fig_pTTPL_beta_CTV_adjusted_association",
    7.5,
    4.5
  )
}

# ------------------------------------------------------------------------------
# 19. SAVE COMPLETE R OBJECT BUNDLE
# ------------------------------------------------------------------------------

saveRDS(
  list(
    input_file =
      input_used,
    design_qc =
      design_qc,
    descriptive =
      descriptive,
    type3 =
      type3_out,
    coefficients =
      coef_out,
    diagnostics =
      diag_out,
    blockfixed =
      fixed_out,
    categorical_vs_continuous =
      model_comparison_out,
    primary_contrasts =
      primary_contrasts,
    primary_overall =
      primary_overall,
    precipitation_slopes =
      precip_slopes,
    all_vs_ctrl =
      vsctrl,
    cohen_effect_specific =
      cohen_out,
    cohen_vs_ctrl =
      cohen_vsctrl,
    cohen_overall =
      overall_cohen,
    fig2a_overall =
      fig2_overall,
    blockmatched =
      blockmatched,
    SEprop_diagnostics =
      meta_diag,
    SEprop_coefficients =
      meta_coefficients,
    SEprop_contrasts =
      meta_contrasts,
    SEprop_overall =
      meta_overall,
    covary_data =
      covary_data,
    covary_results =
      covary_results,
    beta_p_CTV_association =
      ctv_assoc
  ),
  file.path(
    OUTDIR,
    "pTTPL_categorical_primary_overall_Fig2_statistics.rds"
  ),
  version = 3
)

# ------------------------------------------------------------------------------
# 20. README / INTERPRETATION GUIDE
# ------------------------------------------------------------------------------

readme <- c(
  "NWS pTTPL categorical-primary revision statistics",
  "",
  "PRIMARY INPUT",
  input_used,
  "",
  "EXPERIMENTAL STRUCTURE",
  "48 subplot-level trajectories = 12 treatments x 4 blocks per microbial group.",
  "48 subplots are nested within 24 physical whole plots distributed across 4 blocks.",
  "Warming and precipitation are whole-plot factors; clipping is the within-whole-plot split factor.",
  "",
  "PRIMARY MODEL",
  "metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)",
  "WarmF: Unwarmed/Warmed; PrecipF: Ambient/Half/Double; ClipF: Unclipped/Clipped.",
  "Half and Double precipitation are therefore estimated independently rather than constrained to one linear precipitation slope.",
  "",
  "CONTINUOUS-PRECIPITATION SENSITIVITY",
  "metric ~ Warming * Precipitation * Clipping + (1 | Block/WholePlotID)",
  "The numeric model uses Precipitation = -0.5/0/+1 and is retained only as a sensitivity analysis.",
  "03_pTTPL_MODEL_COMPARISON_categorical_vs_continuous.csv compares the two fixed-effect structures using ML fits.",
  "05_pTTPL_SENSITIVITY_continuous_precipitation_slopes.csv reports the continuous-gradient slopes.",
  "",
  "OVERALL FACTOR EFFECTS",
  "05_pTTPL_PRIMARY_categorical_overall_effects_Holm.csv contains four planned overall effects per Group x Metric:",
  "Overall_Warming = equal-weight mean of W-Ctrl, WH-H, WD-D, WC-C, WHC-HC and WDC-DC.",
  "Overall_Half = equal-weight mean of H-Ctrl, WH-W, HC-C and WHC-WC.",
  "Overall_Double = equal-weight mean of D-Ctrl, WD-W, DC-C and WDC-WC.",
  "Overall_Clipping = equal-weight mean of C-Ctrl, WC-W, HC-H, DC-D, WHC-WH and WDC-WD.",
  "P_Holm_four_overall adjusts these four overall P values within each Group x Metric.",
  "",
  "FIG. 2a",
  "07_pTTPL_OVERALL_Cohen_d_Hedges_g.csv contains descriptive overall standardized effect sizes.",
  "08_pTTPL_FIG2a_OVERALL_W_H_D_C_effectsize_and_LMM.csv joins those effect sizes to the corresponding categorical-LMM overall contrasts.",
  "08_pTTPL_FIG2a_OVERALL_beta_compact.csv is the compact beta-only version.",
  "Fig. 2a should use overall Warming/Half/Double/Clipping effects, not W/H/D/C-vs-Ctrl P values.",
  "",
  "TABLE 1",
  "Table 1 remains treatment-combination vs Ctrl.",
  "07_pTTPL_Cohen_d_Hedges_g_vs_Ctrl.csv contains the descriptive treatment-vs-Ctrl Cohen's d / Hedges' g.",
  "06_pTTPL_PRIMARY_all_treatments_vs_Ctrl.csv contains the corresponding categorical-LMM treatment-vs-Ctrl inference.",
  "Do not substitute Fig. 2a overall P values into Table 1 rows.",
  "",
  "MULTIPLE TESTING",
  "Context-specific planned contrasts use Holm adjustment within each Group x Metric x Family.",
  "A global Holm adjustment across all 20 planned contrasts is also supplied.",
  "The four Fig. 2a overall factor effects use Holm adjustment across Warming/Half/Double/Clipping within each Group x Metric.",
  "The 11 treatment-vs-Ctrl contrasts use Holm adjustment across the 11 treatment comparisons within each Group x Metric.",
  "",
  "MODEL DIAGNOSTICS",
  "03_pTTPL_PRIMARY_categorical_model_diagnostics.csv reports random-effect variance components, singularity, marginal/conditional R2 and residual diagnostics.",
  "03_pTTPL_BLOCK_FIXED_categorical_sensitivity_TypeIII.csv treats Block as fixed while retaining the whole-plot random intercept.",
  "",
  "EFFECT SIZE",
  "Cohen's d and Hedges' g are descriptive only.",
  "Formal statistical inference is based on categorical hierarchical LMM contrasts.",
  "Block-matched dz is a conservative sensitivity summary, not the primary inferential test.",
  "",
  "BETA_p SE PROPAGATION",
  "10_pTTPL_Beta_SE_propagation_* files propagate the fitted pTTPL slope sampling variance through V=SE_Beta^2.",
  "The multilevel meta-regression uses WarmF * PrecipF * ClipF, matching the categorical primary treatment structure.",
  "Random effects include Block, whole plot within Block, and plot-level heterogeneity.",
  "CR2 is intentionally not used because there are only four top-level Blocks.",
  "",
  "SLOPE / REFERENCE-VARIANCE CO-VARIATION",
  "12_pTTPL_Beta_LogVar_covary_adjusted_tests.csv uses:",
  "Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID).",
  "H0: gamma = coefficient(Delta_Beta) = 0.",
  "",
  "OPTIONAL COMMUNITY VARIABILITY ASSOCIATION",
  "13_pTTPL_beta_CTV_association_adjusted_descriptive.csv is generated only when the existing",
  "CTV_Hellinger_RMS file is available. It uses:",
  "CTV_Hellinger_RMS ~ Beta_p_raw_z + WarmF * PrecipF * ClipF + (1 | Block/WholePlotID).",
  "",
  "MOST IMPORTANT FILES FOR MANUSCRIPT / RESPONSE",
  "03_pTTPL_PRIMARY_categorical_hierarchical_TypeIII.csv",
  "03_pTTPL_PRIMARY_categorical_model_diagnostics.csv",
  "03_pTTPL_MODEL_COMPARISON_categorical_vs_continuous.csv",
  "04_pTTPL_PRIMARY_categorical_effect_specific_contrasts_Holm.csv",
  "05_pTTPL_PRIMARY_categorical_overall_effects_Holm.csv",
  "06_pTTPL_PRIMARY_all_treatments_vs_Ctrl.csv",
  "07_pTTPL_Cohen_d_Hedges_g_vs_Ctrl.csv",
  "07_pTTPL_OVERALL_Cohen_d_Hedges_g.csv",
  "08_pTTPL_FIG2a_OVERALL_beta_compact.csv",
  "10_pTTPL_Beta_SE_propagation_overall_directional_effects.csv",
  "12_pTTPL_Beta_LogVar_covary_adjusted_tests.csv",
  "14_pTTPL_REVIEWER_Fig2a_OVERALL_effects.csv",
  "14_pTTPL_REVIEWER_Table1_Fig2_vsCtrl_inference.csv",
  "",
  "Do not use a rounded Supplementary Table 2 as statistical input.",
  "The full-precision refit is the authoritative beta_p source."
)

writeLines(
  readme,
  file.path(
    OUTDIR,
    "README_pTTPL_statistical_interpretation.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTDIR,
    "sessionInfo_pTTPL_statistics.txt"
  )
)

message(
  "\n============================================================"
)
message(
  "pTTPL categorical-primary downstream statistics finished."
)
message(
  "Output directory:\n  ",
  normalizePath(
    OUTDIR,
    mustWork = FALSE
  )
)
message(
  "\nFirst inspect:"
)
message(
  "  00_pTTPL_design_QC.csv"
)
message(
  "  03_pTTPL_PRIMARY_categorical_model_diagnostics.csv"
)
message(
  "  03_pTTPL_MODEL_COMPARISON_categorical_vs_continuous.csv"
)
message(
  "  05_pTTPL_PRIMARY_categorical_overall_effects_Holm.csv"
)
message(
  "  08_pTTPL_FIG2a_OVERALL_beta_compact.csv"
)
message(
  "  06_pTTPL_PRIMARY_all_treatments_vs_Ctrl.csv"
)
message(
  "  10_pTTPL_Beta_SE_propagation_overall_directional_effects.csv"
)
message(
  "============================================================\n"
)
