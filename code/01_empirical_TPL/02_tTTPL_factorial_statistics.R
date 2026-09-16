# ==============================================================================
# Reviewer 2 REVISED master analysis (v3)
# Categorical-precipitation primary model + continuous-gradient sensitivity
# Overall-factor Fig. 2 effect sizes + hierarchical inference
# Effect-specific statistics for the 12-treatment x 4-block split-plot design
#
# Experimental/statistical structure used here
# --------------------------------------------
# For each microbial group:
#   48 subplot-level observations
#   = 12 treatment combinations x 4 blocks
#
# Treatment codes:
#   Ctrl W H D C WH WD WC HC DC WHC WDC
#
# Treatment coding used in the PRIMARY model:
#   WarmF   : Unwarmed / Warmed
#   PrecipF : Ambient / Half / Double (categorical; PRIMARY)
#   ClipF   : Unclipped / Clipped
#
# Numeric variables are retained for the continuous-gradient sensitivity:
#   Warming       : 0 = unwarmed, 1 = warmed
#   Precipitation : -0.5 = half, 0 = ambient, 1 = double
#   Clipping      : 0 = unclipped, 1 = clipped
#
# The categorical primary model estimates Half and Double precipitation
# independently. The numeric precipitation model is retained only to test
# whether the simpler linear-gradient restriction is adequate.
#
# The focal simple-effect families are defined by changing ONLY the focal
# treatment factor while holding the other factors constant:
#
# WARMING:
#   W-Ctrl, WH-H, WD-D, WC-C, WHC-HC, WDC-DC
#
# CLIPPING:
#   C-Ctrl, WC-W, HC-H, DC-D, WHC-WH, WDC-WD
#
# HALF / reduced precipitation:
#   H-Ctrl, WH-W, HC-C, WHC-WC
#
# DOUBLE / increased precipitation:
#   D-Ctrl, WD-W, DC-C, WDC-WC
#
# Statistical hierarchy
# ----------------------
# 1. PRIMARY inference:
#      metric ~ WarmF * PrecipF * ClipF +
#               (1 | Block/WholePlotID)
#
#    Precipitation is a three-level categorical experimental factor so that
#    Half-vs-Ambient and Double-vs-Ambient effects are estimated independently.
#
# 2. CONTINUOUS-GRADIENT sensitivity:
#      metric ~ Warming * Precipitation * Clipping +
#               (1 | Block/WholePlotID)
#
#    This model imposes a linear response across precipitation multipliers
#    0.5, 1, and 2 (centered as -0.5, 0, 1). It is compared with the
#    categorical primary model using ML likelihood-ratio tests, AIC/BIC,
#    Nakagawa-style marginal/conditional R2, diagnostics, and optional
#    leave-one-block-out prediction.
#
# 3. Planned effect-specific contrasts:
#    All biological contrasts are obtained from the SAME categorical primary
#    mixed model. Holm correction is reported within each effect family and across
#    all 20 context-specific contrasts. A secondary W/H/D/C-vs-Ctrl table is also
#    retained for diagnostic purposes, but it is NOT the primary Fig. 2 inference.
#
# 4. Fig. 2a OVERALL factor effects:
#    The primary Fig. 2a estimands are overall Warming, Half precipitation,
#    Double precipitation, and Clipping effects, averaged equally across the
#    relevant factorial contexts. For each overall factor effect:
#      - descriptive magnitude = marginal pooled-SD Cohen's d / Hedges' g
#      - formal inference = categorical LMM overall contrast
#      - multiplicity = Holm correction across the four overall effects within
#        each Group x Metric.
#
# 5. Effect-specific magnitude:
#      conventional context-specific Cohen's d / Hedges' g
#      + block-matched standardized difference (dz) as a sensitivity summary.
#
# 6. Legacy response-ratio framework:
#      ordinary response ratio and block-matched log response ratio.
#      Response ratios are calculated only on positive raw-scale quantities.
#      For LogVar_at_ref_raw, RR is calculated on 10^LogVar_at_ref_raw.
#
# 7. Beta-estimation uncertainty:
#      rma.mv with V = SE_Beta^2 and the same categorical fixed-effect structure
#      (WarmF * PrecipF * ClipF), with random intercepts for Block, whole plot
#      nested within Block, and plot-level heterogeneity. Inference uses
#      model-based t tests with containment degrees of freedom. This is a
#      SENSITIVITY ANALYSIS, not a replacement for the primary model.
#
# 8. Figures:
#      - 12-treatment Beta_raw overview
#      - raw vs detrended beta
#      - overall warming effects on raw/detrended beta
#      - warming effect on reference-abundance variance
#      - warming effect on Hellinger community temporal variability
#      - Fig. 2a overall-factor Cohen's d with LMM Holm significance
#      - optional Fig. 2a overall LMM contrast plot
#      - effect-specific forest plots for Beta_raw
#
# 9. Reviewer-specific Beta/variance co-variation analysis:
#      - block-matched warming-induced changes in Beta_raw and LogVar_at_ref_raw
#      - adjusted within-group test of whether Delta_Beta predicts Delta_LogVar_ref
#      - descriptive Pearson/Spearman associations and QC
#      - supplementary delta-delta scatterplot
#
# IMPORTANT
# ---------
# This script re-runs the DOWNSTREAM treatment statistics and reviewer figures.
# It does NOT need to recompute upstream ecological metrics such as Beta_raw,
# Beta_detrended, LogVar_at_ref_raw or CTV_Hellinger_RMS if the current
# plot-level metric table is already correct.
#
# ==============================================================================

rm(list = ls())

# -----------------------------
# 0. User settings
# -----------------------------

# ------------------------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------------------------
# Run this script from the repository root, or set TPL_PROJECT_ROOT to the
# project directory that contains the data folders documented in README.md.
PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

INPUT_FILE <- file.path(
  "results", "01_tTTPL",
  "UPARSE_all_groups_plot_level_revision_metrics.csv"
)

# Fallback if the file was copied to the project root.


OUTDIR <- "./results/02_tTTPL_factorial_statistics"
FIGDIR <- file.path(OUTDIR, "figures")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# IMPORTANT: overwrite policy
# ------------------------------------------------------------------------------
# Keep upstream ecological metric files untouched.
# Remove ONLY obsolete downstream Reviewer-2 statistical outputs generated by
# earlier trial/debug versions, so that the directory contains one final set.
#
# If an old file does not exist, file.remove() simply has nothing to do.
obsolete_files <- c(
  "01_blockfixed_typeIII_ANOVA.csv",
  "06_simple_effects_compact.csv",
  "07_warming_marginal_contrasts_all_metrics.csv",
  "08_hierarchical_factorial_TypeIII_ANOVA.csv",
  "09_beta_variability_association_models.csv",
  "10_hierarchical_context_specific_warming_Holm.csv",
  "12_hierarchical_model_diagnostics.csv",
  "13_fixed_only_sensitivity_TypeIII.csv",
  "14_Beta_SE_propagation_all_treatments_vs_Ctrl.csv",
  "14_Beta_SE_propagation_context_specific_warming.csv",
  "14_Beta_SE_propagation_meta_coefficients.csv",
  "15_UNIFIED_all_treatments_vs_Ctrl.csv",
  "16_UNIFIED_context_specific_warming.csv",
  "17_overall_warming_compact_for_manuscript.csv",
  "18_Beta_SEprop_corrected_model_diagnostics.csv",
  "18_Beta_SEprop_corrected_coefficients.csv",
  "18_Beta_SEprop_corrected_all_treatments_vs_Ctrl.csv",
  "18_Beta_SEprop_corrected_context_specific_warming.csv",
  "18_Beta_SEprop_corrected_overall_warming.csv",
  "18_Beta_SEprop_corrected_interaction_tests.csv",
  "19_Beta_SE_CR2_model_diagnostics.csv",
  "19_Beta_SE_CR2_coefficients.csv",
  "19_Beta_SE_CR2_overall_warming.csv",
  "19_Beta_SE_CR2_context_specific_warming.csv",
  "19_Beta_SE_CR2_all_treatments_vs_Ctrl.csv",
  "19_Beta_SE_CR2_interaction_tests.csv",
  "19_Beta_SE_CR2_comparison_with_nonrobust.csv"
)

obsolete_paths <- file.path(OUTDIR, obsolete_files)
existing_obsolete <- obsolete_paths[file.exists(obsolete_paths)]

if (length(existing_obsolete) > 0) {
  file.remove(existing_obsolete)
}

cat(
  "Removed ", length(existing_obsolete),
  " obsolete Reviewer-2 downstream output file(s).\n",
  sep = ""
)

# Optional reproduction of the original ieggr::compv workflow.
RUN_LEGACY_IEGGR <- TRUE

GROUP_ORDER <- c(
  "Archaea", "Bacteria", "Fungi", "Protists", "Nematoda"
)

TREATMENT_LEVELS <- c(
  "Ctrl", "W", "H", "D", "C",
  "WH", "WD", "WC", "HC", "DC", "WHC", "WDC"
)

# Reviewer-focused metrics already generated by the upstream revision pipeline.
# The script automatically skips a metric if the input table does not contain it.
MAIN_METRICS <- c(
  "Beta_raw",
  "Beta_detrended",
  "LogVar_at_ref_raw",
  "CTV_Hellinger_RMS",
  "WeightedTaxonCV_raw",
  "MedianTaxonCV_raw",
  "MeanPairwise_Hellinger",
  "MeanAnnual_Hellinger",
  "MeanPairwise_BrayCurtis",
  "MeanAnnual_BrayCurtis"
)

BETA_SE_PAIRS <- c(
  Beta_raw       = "SE_Beta_raw",
  Beta_detrended = "SE_Beta_detrended"
)

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "lme4",
  "lmerTest",
  "emmeans",
  "metafor"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Please install required packages first:\n  ",
    paste(missing_pkgs, collapse = ", "),
    "\n\nExample:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(metafor)
})

# Sum contrasts are used globally. Primary treatment inference is obtained from
# explicitly defined estimated marginal means/contrasts, so interpretation does
# not depend on reading lower-order coefficient labels from the raw parameterization.
old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))

emmeans::emm_options(lmer.df = "satterthwaite")

# ------------------------------------------------------------------------------
# 2. CONTRAST DEFINITIONS
# ------------------------------------------------------------------------------

contrast_def <- bind_rows(
  
  # ---------------------------------------------------------------------------
  # Warming effects
  # ---------------------------------------------------------------------------
  tibble(
    Family = "Warming",
    Contrast_order = 1:6,
    Target = c("W", "WH", "WD", "WC", "WHC", "WDC"),
    Reference = c("Ctrl", "H", "D", "C", "HC", "DC"),
    Context = c(
      "Ambient | Unclipped",
      "Half | Unclipped",
      "Double | Unclipped",
      "Ambient | Clipped",
      "Half | Clipped",
      "Double | Clipped"
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Clipping effects
  # Positive estimate = clipped - unclipped
  # ---------------------------------------------------------------------------
  tibble(
    Family = "Clipping",
    Contrast_order = 1:6,
    Target = c("C", "WC", "HC", "DC", "WHC", "WDC"),
    Reference = c("Ctrl", "W", "H", "D", "WH", "WD"),
    Context = c(
      "Unwarmed | Ambient",
      "Warmed | Ambient",
      "Unwarmed | Half",
      "Unwarmed | Double",
      "Warmed | Half",
      "Warmed | Double"
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Reduced precipitation / Half effects
  # Positive estimate = Half - Ambient
  # ---------------------------------------------------------------------------
  tibble(
    Family = "Half_precipitation",
    Contrast_order = 1:4,
    Target = c("H", "WH", "HC", "WHC"),
    Reference = c("Ctrl", "W", "C", "WC"),
    Context = c(
      "Unwarmed | Unclipped",
      "Warmed | Unclipped",
      "Unwarmed | Clipped",
      "Warmed | Clipped"
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Increased precipitation / Double effects
  # Positive estimate = Double - Ambient
  # ---------------------------------------------------------------------------
  tibble(
    Family = "Double_precipitation",
    Contrast_order = 1:4,
    Target = c("D", "WD", "DC", "WDC"),
    Reference = c("Ctrl", "W", "C", "WC"),
    Context = c(
      "Unwarmed | Unclipped",
      "Warmed | Unclipped",
      "Unwarmed | Clipped",
      "Warmed | Clipped"
    )
  )
) %>%
  mutate(
    Contrast = paste0(Target, " - ", Reference),
    Pairing_type = ifelse(
      Family == "Clipping",
      "Within-whole-plot",
      "Block-matched"
    )
  )

write.csv(
  contrast_def,
  file.path(OUTDIR, "00_effect_specific_contrast_definitions.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------------------
# 3. GENERAL HELPERS
# ------------------------------------------------------------------------------

write_out <- function(x, filename) {
  write.csv(
    x,
    file.path(OUTDIR, filename),
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    p <= 0.10  ~ "#",
    TRUE       ~ "ns"
  )
}

make_treatment <- function(warm, precip, clip) {
  
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
  
  unname(map[key])
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

# Equal-weight mean of a collection of simple contrasts.
make_average_weight <- function(simple_weights) {
  Reduce("+", simple_weights) / length(simple_weights)
}

# ------------------------------------------------------------------------------
# 4. READ AND STANDARDIZE DATA
# ------------------------------------------------------------------------------

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found:\n  ", INPUT_FILE,
    "\nPlease edit INPUT_FILE near the top of the script."
  )
}

dat <- read.csv(
  INPUT_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_core <- c(
  "Group", "Plot", "block", "WholePlotID"
)

missing_core <- setdiff(required_core, names(dat))

if (length(missing_core) > 0) {
  stop(
    "Missing core columns: ",
    paste(missing_core, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# NUMERIC treatment coding retained for continuous-gradient sensitivity
# ------------------------------------------------------------------------------
# Warming       : 0 = unwarmed, 1 = warmed
# Precipitation : -0.5 = half, 0 = ambient, 1 = double
# Clipping      : 0 = unclipped, 1 = clipped
#
# Factor variables are used in the PRIMARY categorical model and are also retained
# for treatment labels / descriptive summaries. Numeric variables are retained for
# the continuous-gradient sensitivity model.

# Warming.
if ("WarmF" %in% names(dat)) {
  warm_chr <- as.character(dat$WarmF)
  dat$Warming <- dplyr::case_when(
    warm_chr %in% c("Unwarmed", "Ambient", "0") ~ 0,
    warm_chr %in% c("Warmed", "1") ~ 1,
    TRUE ~ NA_real_
  )
} else if ("warm" %in% names(dat)) {
  warm_num <- safe_num(dat$warm)
  dat$Warming <- dplyr::case_when(
    warm_num == 0 ~ 0,
    warm_num == 1 ~ 1,
    TRUE ~ NA_real_
  )
} else {
  stop("Need either WarmF or warm.")
}

# Precipitation. Accept both the historical raw coding (0.5, 1, 2) and
# the centered coding (-0.5, 0, 1), but STANDARDIZE internally to -0.5/0/1.
if ("PrecipF" %in% names(dat)) {
  precip_chr <- as.character(dat$PrecipF)
  dat$Precipitation <- dplyr::case_when(
    precip_chr %in% c("Half", "half", "-0.5") ~ -0.5,
    precip_chr %in% c("Ambient", "Normal", "normal", "0") ~ 0,
    precip_chr %in% c("Double", "double", "1") ~ 1,
    TRUE ~ NA_real_
  )
} else if ("precip" %in% names(dat)) {
  precip_num <- safe_num(dat$precip)
  observed_precip <- sort(unique(precip_num[is.finite(precip_num)]))
  
  # Historical coding used in the project: 0.5 = half, 1 = ambient, 2 = double.
  if (all(observed_precip %in% c(0.5, 1, 2))) {
    dat$Precipitation <- dplyr::case_when(
      precip_num == 0.5 ~ -0.5,
      precip_num == 1   ~ 0,
      precip_num == 2   ~ 1,
      TRUE ~ NA_real_
    )
    # Centered coding: -0.5 = half, 0 = ambient, 1 = double.
  } else if (all(observed_precip %in% c(-0.5, 0, 1))) {
    dat$Precipitation <- precip_num
  } else {
    stop(
      "Unrecognized precipitation coding. Expected either ",
      "0.5/1/2 or -0.5/0/1."
    )
  }
} else {
  stop("Need either PrecipF or precip.")
}

# Clipping.
if ("ClipF" %in% names(dat)) {
  clip_chr <- as.character(dat$ClipF)
  dat$Clipping <- dplyr::case_when(
    clip_chr %in% c("Unclipped", "0") ~ 0,
    clip_chr %in% c("Clipped", "1") ~ 1,
    TRUE ~ NA_real_
  )
} else if ("clip" %in% names(dat)) {
  clip_num <- safe_num(dat$clip)
  dat$Clipping <- dplyr::case_when(
    clip_num == 0 ~ 0,
    clip_num == 1 ~ 1,
    TRUE ~ NA_real_
  )
} else {
  stop("Need either ClipF or clip.")
}

if (
  any(!is.finite(dat$Warming)) ||
  any(!is.finite(dat$Precipitation)) ||
  any(!is.finite(dat$Clipping))
) {
  bad <- dat %>%
    filter(
      !is.finite(Warming) |
        !is.finite(Precipitation) |
        !is.finite(Clipping)
    )
  
  write_out(
    bad,
    "00_BAD_numeric_treatment_coding.csv"
  )
  
  stop(
    "Some treatment rows could not be standardized to ",
    "Warming=0/1, Precipitation=-0.5/0/1, Clipping=0/1. ",
    "See 00_BAD_numeric_treatment_coding.csv."
  )
}

# Factor variables used by the PRIMARY categorical model and by labels/plots.
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

dat$Group <- factor(
  as.character(dat$Group),
  levels = GROUP_ORDER
)

# Random-effect grouping factors.
dat$Block <- factor(dat$block)
dat$block <- dat$Block
dat$WholePlotID <- factor(dat$WholePlotID)
dat$Plot <- factor(dat$Plot)

# Explicit nested whole-plot identifier. This is useful for diagnostics and
# metafor, even when WholePlotID is already globally unique.
dat$BlockWholePlotID <- interaction(
  dat$Block,
  dat$WholePlotID,
  drop = TRUE,
  sep = ":"
)

# Reconstruct the 12-level treatment factor from the factorial variables.
dat$Treatment <- make_treatment(
  dat$WarmF,
  dat$PrecipF,
  dat$ClipF
)

if (any(is.na(dat$Treatment))) {
  bad <- dat %>%
    filter(is.na(Treatment)) %>%
    select(
      Group, Plot, block, WholePlotID,
      WarmF, PrecipF, ClipF
    )
  
  write_out(
    bad,
    "00_BAD_treatment_reconstruction.csv"
  )
  
  stop(
    "Treatment reconstruction failed for some rows. ",
    "See 00_BAD_treatment_reconstruction.csv."
  )
}

dat$Treatment <- factor(
  dat$Treatment,
  levels = TREATMENT_LEVELS
)

# Keep only available reviewer metrics.
MAIN_METRICS <- intersect(MAIN_METRICS, names(dat))

if (length(MAIN_METRICS) == 0) {
  stop("None of the expected reviewer metrics were found in the input.")
}

# Derived raw-scale reference variance for response-ratio analyses.
if ("LogVar_at_ref_raw" %in% names(dat)) {
  dat$Var_at_ref_raw <- 10^(dat$LogVar_at_ref_raw)
}

# ------------------------------------------------------------------------------
# 5. DESIGN QUALITY CONTROL
# ------------------------------------------------------------------------------

design_qc <- dat %>%
  distinct(
    Group, Plot, Block, WholePlotID, BlockWholePlotID,
    Warming, Precipitation, Clipping,
    WarmF, PrecipF, ClipF, Treatment
  ) %>%
  count(Group, Block, Treatment, name = "N") %>%
  rename(block = Block)

write_out(
  design_qc,
  "01_design_QC_12treatments_x_4blocks.csv"
)

design_qc_summary <- dat %>%
  distinct(
    Group, Plot, Block, BlockWholePlotID, Treatment
  ) %>%
  group_by(Group) %>%
  summarise(
    n_rows = n(),
    n_blocks = n_distinct(Block),
    n_nested_wholeplots = n_distinct(BlockWholePlotID),
    n_treatments = n_distinct(Treatment),
    .groups = "drop"
  ) %>%
  left_join(
    design_qc %>%
      group_by(Group) %>%
      summarise(
        min_n_per_block_treatment = min(N),
        max_n_per_block_treatment = max(N),
        .groups = "drop"
      ),
    by = "Group"
  )

write_out(
  design_qc_summary,
  "01_design_QC_summary.csv"
)

# Validate the random-effects hierarchy used in the STANDARD model.
nested_qc <- dat %>%
  distinct(
    Group,
    Block,
    WholePlotID,
    BlockWholePlotID,
    Warming,
    Precipitation,
    Clipping,
    Treatment
  ) %>%
  group_by(Group, BlockWholePlotID) %>%
  summarise(
    N_subplots = n(),
    N_blocks = n_distinct(Block),
    N_wholeplot_labels = n_distinct(WholePlotID),
    N_clipping_levels = n_distinct(Clipping),
    .groups = "drop"
  ) %>%
  mutate(
    Valid_nested_wholeplot =
      N_blocks == 1 &
      N_wholeplot_labels == 1 &
      N_clipping_levels == 2
  )

write_out(
  nested_qc,
  "01_nested_random_structure_QC.csv"
)

# Whole-plot pairing validation for clipping contrasts.
clipping_qc <- lapply(
  unique(as.character(dat$Group)),
  function(g) {
    
    dg <- dat %>% filter(Group == g)
    
    bind_rows(lapply(
      which(contrast_def$Family == "Clipping"),
      function(j) {
        
        cd <- contrast_def[j, ]
        
        a <- dg %>%
          filter(Treatment == cd$Target) %>%
          transmute(
            block,
            Target = as.character(Treatment),
            WholePlot_target = as.character(WholePlotID)
          )
        
        b <- dg %>%
          filter(Treatment == cd$Reference) %>%
          transmute(
            block,
            Reference = as.character(Treatment),
            WholePlot_reference = as.character(WholePlotID)
          )
        
        z <- inner_join(a, b, by = "block")
        
        tibble(
          Group = g,
          Contrast = cd$Contrast,
          n_blocks = nrow(z),
          all_same_wholeplot =
            nrow(z) > 0 &&
            all(z$WholePlot_target == z$WholePlot_reference)
        )
      }
    ))
  }
) %>%
  bind_rows()

write_out(
  clipping_qc,
  "01_clipping_pairing_QC.csv"
)

# ------------------------------------------------------------------------------
# 6. DESCRIPTIVE TREATMENT SUMMARIES
# ------------------------------------------------------------------------------

descriptive <- bind_rows(lapply(
  MAIN_METRICS,
  function(metric) {
    
    dat %>%
      filter(is.finite(.data[[metric]])) %>%
      group_by(Group, Treatment) %>%
      summarise(
        Metric = metric,
        N = n(),
        Mean = mean(.data[[metric]], na.rm = TRUE),
        SD = sd(.data[[metric]], na.rm = TRUE),
        SE = SD / sqrt(N),
        Median = median(.data[[metric]], na.rm = TRUE),
        .groups = "drop"
      )
  }
))

write_out(
  descriptive,
  "02_descriptive_12treatments_all_metrics.csv"
)

# ------------------------------------------------------------------------------
# 7. PRIMARY CATEGORICAL MODEL + CONTINUOUS-GRADIENT SENSITIVITY
# ------------------------------------------------------------------------------

PRIMARY_MODEL_LABEL <-
  "metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)"

CONTINUOUS_MODEL_LABEL <-
  "metric ~ Warming * Precipitation * Clipping + (1 | Block/WholePlotID)"

# Metrics for which diagnostic plots are saved. Model-comparison tables are
# produced for every available MAIN_METRIC.
DIAGNOSTIC_METRICS <- intersect(
  c("Beta_raw", "Beta_detrended"),
  MAIN_METRICS
)

RUN_LOBO_MODEL_COMPARISON <- TRUE

fit_primary_models <- function(d, metric) {
  
  dd <- d %>%
    filter(is.finite(.data[[metric]])) %>%
    droplevels()
  
  if (nrow(dd) < 20) {
    return(list(
      categorical = NULL,
      continuous = NULL,
      categorical_ml = NULL,
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
  
  # REML fits are used for primary parameter estimation/inference.
  m_cat <- tryCatch(
    lmerTest::lmer(
      form_cat,
      data = dd,
      REML = TRUE
    ),
    error = function(e) e
  )
  
  m_num <- tryCatch(
    lmerTest::lmer(
      form_num,
      data = dd,
      REML = TRUE
    ),
    error = function(e) e
  )
  
  # ML fits are used ONLY to compare the alternative fixed-effect structures.
  m_cat_ml <- tryCatch(
    lme4::lmer(
      form_cat,
      data = dd,
      REML = FALSE
    ),
    error = function(e) e
  )
  
  m_num_ml <- tryCatch(
    lme4::lmer(
      form_num,
      data = dd,
      REML = FALSE
    ),
    error = function(e) e
  )
  
  err <- c()
  if (inherits(m_cat, "error")) {
    err <- c(err, paste0("categorical REML: ", conditionMessage(m_cat)))
    m_cat <- NULL
  }
  if (inherits(m_num, "error")) {
    err <- c(err, paste0("continuous REML: ", conditionMessage(m_num)))
    m_num <- NULL
  }
  if (inherits(m_cat_ml, "error")) {
    err <- c(err, paste0("categorical ML: ", conditionMessage(m_cat_ml)))
    m_cat_ml <- NULL
  }
  if (inherits(m_num_ml, "error")) {
    err <- c(err, paste0("continuous ML: ", conditionMessage(m_num_ml)))
    m_num_ml <- NULL
  }
  
  list(
    categorical = m_cat,
    continuous = m_num,
    categorical_ml = m_cat_ml,
    continuous_ml = m_num_ml,
    data = dd,
    error = if (length(err) == 0) NA_character_ else paste(err, collapse = " | ")
  )
}

extract_type3 <- function(model, group, metric, model_label) {
  an <- as.data.frame(
    stats::anova(
      model,
      type = "III",
      ddf = "Satterthwaite"
    )
  )
  an$Term <- rownames(an)
  rownames(an) <- NULL
  an %>%
    mutate(
      Group = group,
      Metric = metric,
      Model = model_label,
      .before = 1
    )
}

extract_fixed_coefficients <- function(model, group, metric, model_label) {
  co <- as.data.frame(summary(model)$coefficients)
  co$Term <- rownames(co)
  rownames(co) <- NULL
  co %>%
    mutate(
      Group = group,
      Metric = metric,
      Model = model_label,
      .before = 1
    )
}

# Nakagawa-style R2 for Gaussian random-intercept LMMs.
# For the present models, random-effect variance is the sum of random-intercept
# variances; fixed-effect variance is Var(X beta); residual variance is sigma^2.
nakagawa_r2_gaussian <- function(model) {
  if (is.null(model)) {
    return(c(
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      Fixed_variance = NA_real_,
      Random_variance = NA_real_,
      Residual_variance = NA_real_
    ))
  }
  
  fixed_pred <- tryCatch(
    as.numeric(stats::predict(model, re.form = NA)),
    error = function(e) rep(NA_real_, stats::nobs(model))
  )
  fixed_var <- if (sum(is.finite(fixed_pred)) > 1) {
    stats::var(fixed_pred, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  vc <- as.data.frame(lme4::VarCorr(model))
  rand_var <- sum(
    vc$vcov[vc$grp != "Residual"],
    na.rm = TRUE
  )
  resid_var <- vc %>%
    filter(grp == "Residual") %>%
    pull(vcov)
  if (length(resid_var) == 0) resid_var <- NA_real_
  resid_var <- resid_var[1]
  
  total <- fixed_var + rand_var + resid_var
  
  c(
    R2_marginal = ifelse(is.finite(total) && total > 0,
                         fixed_var / total, NA_real_),
    R2_conditional = ifelse(is.finite(total) && total > 0,
                            (fixed_var + rand_var) / total, NA_real_),
    Fixed_variance = fixed_var,
    Random_variance = rand_var,
    Residual_variance = resid_var
  )
}

extract_diag <- function(model, group, metric, n, data_used, model_label) {
  if (is.null(model)) {
    return(tibble(
      Group = group,
      Metric = metric,
      Model = model_label,
      N = n,
      N_blocks = n_distinct(data_used$Block),
      N_nested_wholeplots = n_distinct(data_used$BlockWholePlotID),
      Singular = NA,
      Block_variance = NA_real_,
      WholePlot_within_Block_variance = NA_real_,
      Residual_variance = NA_real_,
      R2_marginal = NA_real_,
      R2_conditional = NA_real_,
      Residual_mean = NA_real_,
      Residual_SD = NA_real_,
      Max_abs_standardized_residual = NA_real_,
      N_abs_standardized_residual_gt3 = NA_integer_,
      QQ_correlation = NA_real_,
      AbsResidual_vs_Fitted_P = NA_real_
    ))
  }
  
  vc <- as.data.frame(lme4::VarCorr(model))
  block_var <- vc %>% filter(grp == "Block") %>% pull(vcov)
  if (length(block_var) == 0) block_var <- NA_real_
  
  wp_var <- vc %>%
    filter(grepl("Block", grp) & grepl("WholePlotID", grp)) %>%
    pull(vcov)
  if (length(wp_var) == 0) wp_var <- NA_real_
  
  resid_var <- vc %>% filter(grp == "Residual") %>% pull(vcov)
  if (length(resid_var) == 0) resid_var <- NA_real_
  
  r2 <- nakagawa_r2_gaussian(model)
  res <- stats::residuals(model)
  fitv <- stats::fitted(model)
  stdres <- res / stats::sigma(model)
  qq_theory <- stats::qnorm(stats::ppoints(length(stdres)))
  qq_obs <- sort(stdres)
  qq_cor <- if (length(stdres) > 2 && stats::sd(stdres) > 0) {
    stats::cor(qq_theory, qq_obs)
  } else {
    NA_real_
  }
  
  het_p <- tryCatch({
    summary(stats::lm(abs(stdres) ~ fitv))$coefficients[2, 4]
  }, error = function(e) NA_real_)
  
  tibble(
    Group = group,
    Metric = metric,
    Model = model_label,
    N = n,
    N_blocks = n_distinct(data_used$Block),
    N_nested_wholeplots = n_distinct(data_used$BlockWholePlotID),
    Singular = lme4::isSingular(model, tol = 1e-6),
    Block_variance = block_var[1],
    WholePlot_within_Block_variance = wp_var[1],
    Residual_variance = resid_var[1],
    R2_marginal = unname(r2["R2_marginal"]),
    R2_conditional = unname(r2["R2_conditional"]),
    Residual_mean = mean(res, na.rm = TRUE),
    Residual_SD = sd(res, na.rm = TRUE),
    Max_abs_standardized_residual = max(abs(stdres), na.rm = TRUE),
    N_abs_standardized_residual_gt3 = sum(abs(stdres) > 3, na.rm = TRUE),
    QQ_correlation = qq_cor,
    AbsResidual_vs_Fitted_P = het_p
  )
}

# Compare the categorical primary model with the constrained continuous-gradient
# model. Both models are fitted by ML here because their fixed-effect structures
# differ. The categorical model contains the continuous-gradient model as a
# constrained submodel for the three precipitation levels.
compare_cat_vs_continuous <- function(m_cat_ml, m_num_ml, group, metric) {
  if (is.null(m_cat_ml) || is.null(m_num_ml)) {
    return(tibble(
      Group = group,
      Metric = metric,
      N = NA_integer_,
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
      LR_P = NA_real_,
      Preferred_by_AIC = NA_character_
    ))
  }
  
  cmp <- tryCatch(
    as.data.frame(stats::anova(m_num_ml, m_cat_ml, refit = FALSE)),
    error = function(e) NULL
  )
  
  aic_cat <- AIC(m_cat_ml)
  aic_num <- AIC(m_num_ml)
  bic_cat <- BIC(m_cat_ml)
  bic_num <- BIC(m_num_ml)
  ll_cat <- as.numeric(logLik(m_cat_ml))
  ll_num <- as.numeric(logLik(m_num_ml))
  
  lr_chisq <- lr_df <- lr_p <- NA_real_
  if (!is.null(cmp) && nrow(cmp) >= 2) {
    if ("Chisq" %in% names(cmp)) lr_chisq <- cmp$Chisq[2]
    if ("Df" %in% names(cmp)) lr_df <- cmp$Df[2]
    pcol <- grep("Pr\\(>Chisq\\)", names(cmp), value = TRUE)
    if (length(pcol) == 1) lr_p <- cmp[[pcol]][2]
  }
  
  tibble(
    Group = group,
    Metric = metric,
    N = stats::nobs(m_cat_ml),
    AIC_categorical = aic_cat,
    AIC_continuous = aic_num,
    Delta_AIC_cat_minus_cont = aic_cat - aic_num,
    BIC_categorical = bic_cat,
    BIC_continuous = bic_num,
    Delta_BIC_cat_minus_cont = bic_cat - bic_num,
    LogLik_categorical = ll_cat,
    LogLik_continuous = ll_num,
    LR_Chisq = lr_chisq,
    LR_df = lr_df,
    LR_P = lr_p,
    Preferred_by_AIC = case_when(
      is.na(aic_cat) | is.na(aic_num) ~ NA_character_,
      aic_cat + 2 < aic_num ~ "Categorical",
      aic_num + 2 < aic_cat ~ "Continuous",
      TRUE ~ "Similar (Delta AIC < 2)"
    )
  )
}

# Leave-one-block-out prediction compares fixed-effect generalization to a new
# block. Predictions exclude random effects (re.form = NA), because the held-out
# block has no estimated random intercept in the training data.
lobo_compare_models <- function(d, metric) {
  blocks <- unique(as.character(d$Block))
  rows <- list()
  idx <- 1
  
  for (b in blocks) {
    train <- d %>% filter(as.character(Block) != b) %>% droplevels()
    test  <- d %>% filter(as.character(Block) == b)
    
    form_cat <- as.formula(
      paste0(metric, " ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)")
    )
    form_num <- as.formula(
      paste0(metric, " ~ Warming * Precipitation * Clipping + (1 | Block/WholePlotID)")
    )
    
    mcat <- tryCatch(
      lme4::lmer(form_cat, data = train, REML = FALSE),
      error = function(e) NULL
    )
    mnum <- tryCatch(
      lme4::lmer(form_num, data = train, REML = FALSE),
      error = function(e) NULL
    )
    
    y <- test[[metric]]
    pcat <- if (!is.null(mcat)) tryCatch(
      predict(mcat, newdata = test, re.form = NA, allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, nrow(test))
    ) else rep(NA_real_, nrow(test))
    pnum <- if (!is.null(mnum)) tryCatch(
      predict(mnum, newdata = test, re.form = NA, allow.new.levels = TRUE),
      error = function(e) rep(NA_real_, nrow(test))
    ) else rep(NA_real_, nrow(test))
    
    score <- function(pred) {
      ok <- is.finite(y) & is.finite(pred)
      if (sum(ok) == 0) return(c(RMSE = NA_real_, MAE = NA_real_))
      c(
        RMSE = sqrt(mean((y[ok] - pred[ok])^2)),
        MAE = mean(abs(y[ok] - pred[ok]))
      )
    }
    sc_cat <- score(pcat)
    sc_num <- score(pnum)
    
    rows[[idx]] <- tibble(
      Held_out_block = b,
      N_test = nrow(test),
      RMSE_categorical = unname(sc_cat["RMSE"]),
      RMSE_continuous = unname(sc_num["RMSE"]),
      MAE_categorical = unname(sc_cat["MAE"]),
      MAE_continuous = unname(sc_num["MAE"])
    )
    idx <- idx + 1
  }
  
  bind_rows(rows)
}

# Sensitivity model with Block fixed. This is especially useful because Block
# has only four levels; the whole-plot random intercept is retained.
fit_blockfixed_sensitivity <- function(d, metric) {
  form <- as.formula(
    paste0(
      metric,
      " ~ WarmF * PrecipF * ClipF + ",
      "Block + (1 | BlockWholePlotID)"
    )
  )
  lmerTest::lmer(form, data = d, REML = TRUE)
}

all_type3 <- list()
all_coef <- list()
all_diag <- list()
all_fixed <- list()
all_cont_type3 <- list()
all_cont_coef <- list()
all_cont_diag <- list()
all_model_cmp <- list()
all_lobo <- list()
primary_models <- list()

idx <- 1
idx_lobo <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  dg <- dat %>% filter(Group == g) %>% droplevels()
  
  for (metric in MAIN_METRICS) {
    ans <- fit_primary_models(dg, metric)
    key <- paste(g, metric, sep = "::")
    primary_models[[key]] <- ans
    
    if (!is.null(ans$categorical)) {
      all_type3[[idx]] <- extract_type3(
        ans$categorical, g, metric, PRIMARY_MODEL_LABEL
      )
      all_coef[[idx]] <- extract_fixed_coefficients(
        ans$categorical, g, metric, PRIMARY_MODEL_LABEL
      )
      all_diag[[idx]] <- extract_diag(
        ans$categorical, g, metric, nrow(ans$data), ans$data,
        PRIMARY_MODEL_LABEL
      ) %>% mutate(Error = ans$error)
      
      mf <- tryCatch(
        fit_blockfixed_sensitivity(ans$data, metric),
        error = function(e) e
      )
      if (!inherits(mf, "error")) {
        af2 <- as.data.frame(
          stats::anova(mf, type = "III", ddf = "Satterthwaite")
        )
        af2$Term <- rownames(af2)
        rownames(af2) <- NULL
        all_fixed[[idx]] <- af2 %>%
          mutate(
            Group = g,
            Metric = metric,
            Model =
              "metric ~ WarmF * PrecipF * ClipF + Block + (1 | BlockWholePlotID)",
            Singular = lme4::isSingular(mf, tol = 1e-6),
            .before = 1
          )
      }
    } else {
      all_diag[[idx]] <- extract_diag(
        NULL, g, metric, nrow(ans$data), ans$data, PRIMARY_MODEL_LABEL
      ) %>% mutate(Error = ans$error)
    }
    
    if (!is.null(ans$continuous)) {
      all_cont_type3[[idx]] <- extract_type3(
        ans$continuous, g, metric, CONTINUOUS_MODEL_LABEL
      )
      all_cont_coef[[idx]] <- extract_fixed_coefficients(
        ans$continuous, g, metric, CONTINUOUS_MODEL_LABEL
      )
      all_cont_diag[[idx]] <- extract_diag(
        ans$continuous, g, metric, nrow(ans$data), ans$data,
        CONTINUOUS_MODEL_LABEL
      ) %>% mutate(Error = ans$error)
    }
    
    all_model_cmp[[idx]] <- compare_cat_vs_continuous(
      ans$categorical_ml,
      ans$continuous_ml,
      g,
      metric
    )
    
    if (RUN_LOBO_MODEL_COMPARISON && metric %in% DIAGNOSTIC_METRICS) {
      ll <- lobo_compare_models(ans$data, metric) %>%
        mutate(Group = g, Metric = metric, .before = 1)
      all_lobo[[idx_lobo]] <- ll
      idx_lobo <- idx_lobo + 1
    }
    
    idx <- idx + 1
  }
}

type3_out <- bind_rows(all_type3)
coef_out <- bind_rows(all_coef)
diag_out <- bind_rows(all_diag)
fixed_out <- bind_rows(all_fixed)
continuous_type3_out <- bind_rows(all_cont_type3)
continuous_coef_out <- bind_rows(all_cont_coef)
continuous_diag_out <- bind_rows(all_cont_diag)
model_comparison_out <- bind_rows(all_model_cmp)
lobo_out <- bind_rows(all_lobo)

write_out(type3_out, "03_PRIMARY_categorical_hierarchical_TypeIII.csv")
write_out(coef_out, "03_PRIMARY_categorical_model_coefficients.csv")
write_out(diag_out, "03_PRIMARY_categorical_model_diagnostics.csv")
write_out(fixed_out, "03_BLOCK_FIXED_categorical_sensitivity_TypeIII.csv")
write_out(continuous_type3_out, "03_CONTINUOUS_sensitivity_TypeIII.csv")
write_out(continuous_coef_out, "03_CONTINUOUS_sensitivity_coefficients.csv")
write_out(continuous_diag_out, "03_CONTINUOUS_sensitivity_model_diagnostics.csv")
write_out(model_comparison_out, "03_MODEL_COMPARISON_categorical_vs_continuous.csv")
if (nrow(lobo_out) > 0) {
  write_out(lobo_out, "03_MODEL_COMPARISON_leave_one_block_out.csv")
  lobo_summary <- lobo_out %>%
    group_by(Group, Metric) %>%
    summarise(
      Mean_RMSE_categorical = mean(RMSE_categorical, na.rm = TRUE),
      Mean_RMSE_continuous = mean(RMSE_continuous, na.rm = TRUE),
      Mean_MAE_categorical = mean(MAE_categorical, na.rm = TRUE),
      Mean_MAE_continuous = mean(MAE_continuous, na.rm = TRUE),
      Blocks_cat_lower_RMSE = sum(RMSE_categorical < RMSE_continuous, na.rm = TRUE),
      Blocks_cont_lower_RMSE = sum(RMSE_continuous < RMSE_categorical, na.rm = TRUE),
      .groups = "drop"
    )
  write_out(lobo_summary, "03_MODEL_COMPARISON_leave_one_block_out_summary.csv")
}

# Save residual-vs-fitted and QQ plots for the primary categorical beta models.
for (key in names(primary_models)) {
  ans <- primary_models[[key]]
  if (is.null(ans$categorical)) next
  pieces <- strsplit(key, "::", fixed = TRUE)[[1]]
  g <- pieces[1]
  metric <- pieces[2]
  if (!metric %in% DIAGNOSTIC_METRICS) next
  
  dd_diag <- tibble(
    Fitted = fitted(ans$categorical),
    Residual = residuals(ans$categorical),
    StdResidual = residuals(ans$categorical) / sigma(ans$categorical)
  )
  
  p_res <- ggplot(dd_diag, aes(x = Fitted, y = StdResidual)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2, alpha = 0.75) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.6) +
    labs(
      x = "Fitted value",
      y = "Standardized residual",
      title = paste(g, metric, "primary categorical model")
    ) +
    theme_bw(base_size = 11)
  
  p_qq <- ggplot(dd_diag, aes(sample = StdResidual)) +
    stat_qq(size = 2, alpha = 0.75) +
    stat_qq_line(linewidth = 0.6) +
    labs(
      x = "Theoretical quantile",
      y = "Standardized residual",
      title = paste(g, metric, "primary categorical model")
    ) +
    theme_bw(base_size = 11)
  
  safe_key <- gsub("[^A-Za-z0-9]+", "_", paste(g, metric, sep = "_"))
  ggsave(
    file.path(FIGDIR, paste0("Diag_residual_fitted_", safe_key, ".pdf")),
    p_res, width = 5.5, height = 4.2
  )
  ggsave(
    file.path(FIGDIR, paste0("Diag_QQ_", safe_key, ".pdf")),
    p_qq, width = 5.5, height = 4.2
  )
}

# ------------------------------------------------------------------------------
# 8. PLANNED CONTRASTS FROM THE PRIMARY CATEGORICAL MODEL
# ------------------------------------------------------------------------------

make_emm_weight <- function(emm_levels, target, reference) {
  w <- rep(0, length(emm_levels))
  names(w) <- emm_levels
  if (!all(c(target, reference) %in% emm_levels)) {
    stop(
      "Missing treatment level in EMM grid: ",
      paste(setdiff(c(target, reference), emm_levels), collapse = ", ")
    )
  }
  w[target] <- 1
  w[reference] <- -1
  w
}

label_factor_emm_treatment <- function(warmf, precipf, clipf) {
  make_treatment(
    as.character(warmf),
    as.character(precipf),
    as.character(clipf)
  )
}

extract_categorical_contrasts <- function(model, group, metric) {
  emm <- emmeans::emmeans(
    model,
    specs = ~ WarmF * PrecipF * ClipF,
    lmer.df = "satterthwaite"
  )
  emm_grid <- as.data.frame(emm)
  emm_levels <- label_factor_emm_treatment(
    emm_grid$WarmF,
    emm_grid$PrecipF,
    emm_grid$ClipF
  )
  
  if (any(is.na(emm_levels)) || anyDuplicated(emm_levels) > 0 ||
      length(emm_levels) != 12) {
    stop("Could not reconstruct the 12-treatment categorical EMM grid.")
  }
  
  rows <- list()
  weights_by_family <- list()
  
  for (fam in unique(contrast_def$Family)) {
    defs <- contrast_def %>%
      filter(Family == fam) %>%
      arrange(Contrast_order)
    fam_weights <- list()
    
    for (j in seq_len(nrow(defs))) {
      w <- make_emm_weight(emm_levels, defs$Target[j], defs$Reference[j])
      fam_weights[[j]] <- w
      cc <- emmeans::contrast(
        emm,
        method = setNames(list(w), defs$Contrast[j]),
        adjust = "none"
      )
      ss <- as.data.frame(summary(cc, infer = c(TRUE, TRUE), adjust = "none"))
      
      rows[[length(rows) + 1]] <- tibble(
        Group = group,
        Metric = metric,
        Family = fam,
        Contrast_order = defs$Contrast_order[j],
        Contrast = defs$Contrast[j],
        Target = defs$Target[j],
        Reference = defs$Reference[j],
        Context = defs$Context[j],
        Pairing_type = defs$Pairing_type[j],
        Model_basis = PRIMARY_MODEL_LABEL,
        Estimate = ss$estimate[1],
        SE = ss$SE[1],
        df = ss$df[1],
        CI_low = ss$lower.CL[1],
        CI_high = ss$upper.CL[1],
        t_ratio = ss$t.ratio[1],
        P_raw = ss$p.value[1]
      )
    }
    weights_by_family[[fam]] <- fam_weights
  }
  
  individual <- bind_rows(rows) %>%
    group_by(Group, Metric, Family) %>%
    mutate(
      P_Holm_within_family = p.adjust(P_raw, method = "holm"),
      Significance_Holm = sig_label(P_Holm_within_family)
    ) %>%
    ungroup() %>%
    group_by(Group, Metric) %>%
    mutate(P_Holm_all20 = p.adjust(P_raw, method = "holm")) %>%
    ungroup()
  
  overall_rows <- list()
  for (fam in names(weights_by_family)) {
    wavg <- make_average_weight(weights_by_family[[fam]])
    label <- paste0("Overall_", fam)
    cc <- emmeans::contrast(
      emm,
      method = setNames(list(wavg), label),
      adjust = "none"
    )
    ss <- as.data.frame(summary(cc, infer = c(TRUE, TRUE), adjust = "none"))
    overall_rows[[length(overall_rows) + 1]] <- tibble(
      Group = group,
      Metric = metric,
      Family = fam,
      Contrast = label,
      Model_basis = PRIMARY_MODEL_LABEL,
      Estimate = ss$estimate[1],
      SE = ss$SE[1],
      df = ss$df[1],
      CI_low = ss$lower.CL[1],
      CI_high = ss$upper.CL[1],
      t_ratio = ss$t.ratio[1],
      P_raw = ss$p.value[1]
    )
  }
  
  overall <- bind_rows(overall_rows) %>%
    group_by(Group, Metric) %>%
    mutate(
      P_Holm_four_overall = p.adjust(P_raw, method = "holm"),
      Significance_Holm_four_overall = sig_label(P_Holm_four_overall)
    ) %>%
    ungroup()
  
  list(individual = individual, overall = overall)
}

# Continuous sensitivity contrasts are also retained to show how much the
# linear-gradient restriction changes the model-adjusted W/H/D/C estimates.
label_numeric_emm_treatment <- function(warming, precipitation, clipping) {
  wf <- ifelse(warming == 1, "Warmed", "Unwarmed")
  pf <- case_when(
    precipitation == -0.5 ~ "Half",
    precipitation == 0 ~ "Ambient",
    precipitation == 1 ~ "Double",
    TRUE ~ NA_character_
  )
  cf <- ifelse(clipping == 1, "Clipped", "Unclipped")
  make_treatment(wf, pf, cf)
}

extract_continuous_contrasts <- function(model, group, metric) {
  emm <- emmeans::emmeans(
    model,
    specs = ~ Warming * Precipitation * Clipping,
    at = list(
      Warming = c(0, 1),
      Precipitation = c(-0.5, 0, 1),
      Clipping = c(0, 1)
    ),
    lmer.df = "satterthwaite"
  )
  eg <- as.data.frame(emm)
  levels12 <- label_numeric_emm_treatment(
    eg$Warming, eg$Precipitation, eg$Clipping
  )
  rows <- list()
  for (j in seq_len(nrow(contrast_def))) {
    cd <- contrast_def[j, ]
    w <- make_emm_weight(levels12, cd$Target, cd$Reference)
    cc <- emmeans::contrast(
      emm,
      method = setNames(list(w), cd$Contrast),
      adjust = "none"
    )
    ss <- as.data.frame(summary(cc, infer = c(TRUE, TRUE), adjust = "none"))
    rows[[j]] <- tibble(
      Group = group,
      Metric = metric,
      Family = cd$Family,
      Contrast_order = cd$Contrast_order,
      Contrast = cd$Contrast,
      Target = cd$Target,
      Reference = cd$Reference,
      Context = cd$Context,
      Model_basis = CONTINUOUS_MODEL_LABEL,
      Estimate_continuous = ss$estimate[1],
      SE_continuous = ss$SE[1],
      df_continuous = ss$df[1],
      CI_low_continuous = ss$lower.CL[1],
      CI_high_continuous = ss$upper.CL[1],
      P_continuous = ss$p.value[1]
    )
  }
  bind_rows(rows)
}

all_primary_contrasts <- list()
all_primary_overall <- list()
all_continuous_contrasts <- list()
idx <- 1
idx_cont <- 1

for (key in names(primary_models)) {
  ans <- primary_models[[key]]
  pieces <- strsplit(key, "::", fixed = TRUE)[[1]]
  g <- pieces[1]
  metric <- pieces[2]
  
  if (!is.null(ans$categorical)) {
    rr <- tryCatch(
      extract_categorical_contrasts(ans$categorical, g, metric),
      error = function(e) e
    )
    if (!inherits(rr, "error")) {
      all_primary_contrasts[[idx]] <- rr$individual
      all_primary_overall[[idx]] <- rr$overall
      idx <- idx + 1
    } else {
      warning("Categorical contrast extraction failed for ", key, ": ",
              conditionMessage(rr))
    }
  }
  
  if (!is.null(ans$continuous)) {
    cr <- tryCatch(
      extract_continuous_contrasts(ans$continuous, g, metric),
      error = function(e) e
    )
    if (!inherits(cr, "error")) {
      all_continuous_contrasts[[idx_cont]] <- cr
      idx_cont <- idx_cont + 1
    }
  }
}

primary_contrasts <- bind_rows(all_primary_contrasts)
primary_overall <- bind_rows(all_primary_overall)
continuous_contrasts <- bind_rows(all_continuous_contrasts)

# Secondary diagnostic: four single-factor treatment cells compared with Ctrl.
# These are retained because the original Fig. 2 used W/H/D/C-vs-Ctrl Cohen's d.
# They are NOT the primary Fig. 2a estimands in v3; Fig. 2a now uses overall
# factor effects defined below.
fig2_single4 <- primary_contrasts %>%
  filter(Contrast %in% c("W - Ctrl", "H - Ctrl", "D - Ctrl", "C - Ctrl")) %>%
  mutate(
    Fig2_order = match(Contrast, c("W - Ctrl", "H - Ctrl", "D - Ctrl", "C - Ctrl")),
    Fig2_treatment = sub(" - Ctrl$", "", Contrast)
  ) %>%
  group_by(Group, Metric) %>%
  mutate(
    P_Holm_single4 = p.adjust(P_raw, method = "holm"),
    Significance_Holm_single4 = sig_label(P_Holm_single4)
  ) %>%
  ungroup() %>%
  arrange(Metric, Group, Fig2_order)

write_out(primary_contrasts, "04_PRIMARY_effect_specific_contrasts_Holm.csv")
write_out(fig2_single4, "04_PRIMARY_single_treatment_vs_Ctrl_Holm.csv")
write_out(primary_overall, "05_PRIMARY_overall_directional_effects.csv")
write_out(continuous_contrasts, "05_CONTINUOUS_sensitivity_effect_specific_contrasts.csv")

# Direct categorical-vs-continuous contrast comparison.
contrast_model_comparison <- primary_contrasts %>%
  select(
    Group, Metric, Family, Contrast, Target, Reference, Context,
    Estimate_categorical = Estimate,
    SE_categorical = SE,
    P_categorical = P_raw
  ) %>%
  left_join(
    continuous_contrasts %>%
      select(
        Group, Metric, Family, Contrast,
        Estimate_continuous, SE_continuous, P_continuous
      ),
    by = c("Group", "Metric", "Family", "Contrast")
  ) %>%
  mutate(
    Delta_estimate_cat_minus_cont = Estimate_categorical - Estimate_continuous,
    Direction_changed =
      sign(Estimate_categorical) != sign(Estimate_continuous) &
      is.finite(Estimate_categorical) & is.finite(Estimate_continuous)
  )
write_out(
  contrast_model_comparison,
  "05_MODEL_COMPARISON_effect_specific_categorical_vs_continuous.csv"
)

# Continuous precipitation slopes are retained as a sensitivity summary.
extract_precipitation_slopes <- function(model, group, metric) {
  tr <- emmeans::emtrends(
    model,
    specs = ~ Warming * Clipping,
    var = "Precipitation",
    at = list(Warming = c(0, 1), Clipping = c(0, 1)),
    lmer.df = "satterthwaite"
  )
  ss <- as.data.frame(summary(tr, infer = c(TRUE, TRUE), adjust = "none"))
  slope_col <- grep("\\.trend$", names(ss), value = TRUE)
  if (length(slope_col) != 1) {
    stop("Could not identify the precipitation trend column from emtrends().")
  }
  context_rows <- tibble(
    Group = group,
    Metric = metric,
    Scope = "Context_specific",
    Warming = ss$Warming,
    Clipping = ss$Clipping,
    Context = paste0(
      ifelse(ss$Warming == 1, "Warmed", "Unwarmed"),
      " | ",
      ifelse(ss$Clipping == 1, "Clipped", "Unclipped")
    ),
    Precipitation_slope = ss[[slope_col]],
    SE = ss$SE,
    df = ss$df,
    CI_low = ss$lower.CL,
    CI_high = ss$upper.CL,
    t_ratio = ss$t.ratio,
    P = ss$p.value
  )
  w <- rep(1 / nrow(ss), nrow(ss))
  ov <- emmeans::contrast(
    tr,
    method = list(Overall_precipitation_slope = w),
    adjust = "none"
  )
  oo <- as.data.frame(summary(ov, infer = c(TRUE, TRUE), adjust = "none"))
  overall_row <- tibble(
    Group = group,
    Metric = metric,
    Scope = "Overall",
    Warming = NA_real_,
    Clipping = NA_real_,
    Context = "Equal-weight average across warming x clipping contexts",
    Precipitation_slope = oo$estimate[1],
    SE = oo$SE[1],
    df = oo$df[1],
    CI_low = oo$lower.CL[1],
    CI_high = oo$upper.CL[1],
    t_ratio = oo$t.ratio[1],
    P = oo$p.value[1]
  )
  bind_rows(context_rows, overall_row)
}

precip_slope_rows <- list()
idx_ps <- 1
for (key in names(primary_models)) {
  ans <- primary_models[[key]]
  if (is.null(ans$continuous)) next
  pieces <- strsplit(key, "::", fixed = TRUE)[[1]]
  rr <- tryCatch(
    extract_precipitation_slopes(ans$continuous, pieces[1], pieces[2]),
    error = function(e) e
  )
  if (!inherits(rr, "error")) {
    precip_slope_rows[[idx_ps]] <- rr
    idx_ps <- idx_ps + 1
  }
}
precip_slopes <- bind_rows(precip_slope_rows)
if (nrow(precip_slopes) > 0) {
  precip_slopes <- precip_slopes %>%
    group_by(Group, Metric, Scope) %>%
    mutate(P_Holm_within_scope = p.adjust(P, method = "holm")) %>%
    ungroup()
}
write_out(
  precip_slopes,
  "05_CONTINUOUS_sensitivity_precipitation_slopes.csv"
)

# ------------------------------------------------------------------------------
# 9. BLOCK-MATCHED DIFFERENCE SENSITIVITY
#
# For Warming, Half and Double:
#   the two observations are matched by BLOCK but are from different whole plots.
#
# For Clipping:
#   the two observations are both block matched AND originate from the same
#   physical whole plot. This is validated above.
#
# These n=4 block-difference tests are intentionally conservative and are NOT
# used as a replacement for the primary hierarchical model.
# ------------------------------------------------------------------------------

get_matched_pair <- function(d, metric, target, reference) {
  
  a <- d %>%
    filter(Treatment == target) %>%
    transmute(
      block,
      Target_value = .data[[metric]],
      WholePlot_target = as.character(WholePlotID),
      Plot_target = as.character(Plot)
    )
  
  b <- d %>%
    filter(Treatment == reference) %>%
    transmute(
      block,
      Reference_value = .data[[metric]],
      WholePlot_reference = as.character(WholePlotID),
      Plot_reference = as.character(Plot)
    )
  
  inner_join(a, b, by = "block") %>%
    mutate(
      Difference = Target_value - Reference_value,
      SameWholePlot =
        WholePlot_target == WholePlot_reference
    )
}

safe_t_onesample <- function(x) {
  
  x <- x[is.finite(x)]
  
  if (length(x) < 2 || sd(x) == 0) {
    return(list(
      estimate = ifelse(length(x) > 0, mean(x), NA_real_),
      se = NA_real_,
      df = length(x) - 1,
      lo = NA_real_,
      hi = NA_real_,
      p = NA_real_
    ))
  }
  
  tt <- t.test(x, mu = 0)
  
  list(
    estimate = mean(x),
    se = sd(x) / sqrt(length(x)),
    df = unname(tt$parameter),
    lo = tt$conf.int[1],
    hi = tt$conf.int[2],
    p = tt$p.value
  )
}

blockmatched_rows <- list()
idx <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  
  dg <- dat %>% filter(Group == g)
  
  for (metric in MAIN_METRICS) {
    
    for (j in seq_len(nrow(contrast_def))) {
      
      cd <- contrast_def[j, ]
      
      pp <- get_matched_pair(
        dg,
        metric,
        cd$Target,
        cd$Reference
      )
      
      tt <- safe_t_onesample(pp$Difference)
      
      dz <- if (
        nrow(pp) >= 2 &&
        is.finite(sd(pp$Difference)) &&
        sd(pp$Difference) > 0
      ) {
        mean(pp$Difference) / sd(pp$Difference)
      } else {
        NA_real_
      }
      
      blockmatched_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        Family = cd$Family,
        Contrast_order = cd$Contrast_order,
        Contrast = cd$Contrast,
        Target = cd$Target,
        Reference = cd$Reference,
        Context = cd$Context,
        Pairing_type = cd$Pairing_type,
        N_blocks = nrow(pp),
        All_same_wholeplot =
          ifelse(
            cd$Family == "Clipping",
            all(pp$SameWholePlot),
            NA
          ),
        Mean_block_difference = tt$estimate,
        SE_block_difference = tt$se,
        df = tt$df,
        CI_low = tt$lo,
        CI_high = tt$hi,
        P_blockmatched = tt$p,
        Cohen_dz_blockmatched = dz
      )
      
      idx <- idx + 1
    }
  }
}

blockmatched <- bind_rows(blockmatched_rows) %>%
  group_by(Group, Metric, Family) %>%
  mutate(
    P_blockmatched_Holm =
      p.adjust(P_blockmatched, method = "holm"),
    Significance_blockmatched =
      sig_label(P_blockmatched_Holm)
  ) %>%
  ungroup()

write_out(
  blockmatched,
  "06_blockmatched_effect_specific_sensitivity.csv"
)

# Block-matched overall family effect:
# first calculate each simple contrast within each block, then average the
# family-specific contrasts within that block, finally test the four block means.
block_overall_rows <- list()
idx <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  
  dg <- dat %>% filter(Group == g)
  
  for (metric in MAIN_METRICS) {
    
    for (fam in unique(contrast_def$Family)) {
      
      defs <- contrast_def %>%
        filter(Family == fam)
      
      per_contrast <- bind_rows(lapply(
        seq_len(nrow(defs)),
        function(j) {
          
          pp <- get_matched_pair(
            dg,
            metric,
            defs$Target[j],
            defs$Reference[j]
          )
          
          pp %>%
            transmute(
              block,
              Contrast = defs$Contrast[j],
              Difference
            )
        }
      ))
      
      per_block <- per_contrast %>%
        group_by(block) %>%
        summarise(
          Family_mean_difference =
            mean(Difference, na.rm = TRUE),
          .groups = "drop"
        )
      
      tt <- safe_t_onesample(
        per_block$Family_mean_difference
      )
      
      block_overall_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        Family = fam,
        N_blocks = nrow(per_block),
        Estimate = tt$estimate,
        SE = tt$se,
        df = tt$df,
        CI_low = tt$lo,
        CI_high = tt$hi,
        P_blockmatched_overall = tt$p
      )
      
      idx <- idx + 1
    }
  }
}

block_overall <- bind_rows(block_overall_rows) %>%
  group_by(Group, Metric) %>%
  mutate(
    P_Holm_four_overall =
      p.adjust(P_blockmatched_overall, method = "holm")
  ) %>%
  ungroup()

write_out(
  block_overall,
  "06_blockmatched_overall_directional_effects.csv"
)

# ------------------------------------------------------------------------------
# 10. COHEN'S d / HEDGES' g
#
# Conventional pooled-SD d is retained for continuity with the original paper.
# It is descriptive effect magnitude, NOT the primary significance test.
# ------------------------------------------------------------------------------

cohen_pooled <- function(x, y) {
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  n1 <- length(x)
  n0 <- length(y)
  
  if (n1 < 2 || n0 < 2) {
    return(c(
      n_target = n1,
      n_reference = n0,
      mean_target = NA,
      mean_reference = NA,
      sd_target = NA,
      sd_reference = NA,
      sd_pooled = NA,
      d = NA,
      g = NA
    ))
  }
  
  s1 <- sd(x)
  s0 <- sd(y)
  
  sp <- sqrt(
    ((n1 - 1) * s1^2 + (n0 - 1) * s0^2) /
      (n1 + n0 - 2)
  )
  
  d <- if (is.finite(sp) && sp > 0) {
    (mean(x) - mean(y)) / sp
  } else {
    NA_real_
  }
  
  df <- n1 + n0 - 2
  
  # Small-sample Hedges correction.
  J <- if (df > 1) {
    1 - 3 / (4 * df - 1)
  } else {
    NA_real_
  }
  
  g <- J * d
  
  c(
    n_target = n1,
    n_reference = n0,
    mean_target = mean(x),
    mean_reference = mean(y),
    sd_target = s1,
    sd_reference = s0,
    sd_pooled = sp,
    d = d,
    g = g
  )
}

cohen_rows <- list()
idx <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  
  dg <- dat %>% filter(Group == g)
  
  for (metric in MAIN_METRICS) {
    
    for (j in seq_len(nrow(contrast_def))) {
      
      cd <- contrast_def[j, ]
      
      x <- dg %>%
        filter(Treatment == cd$Target) %>%
        pull(all_of(metric))
      
      y <- dg %>%
        filter(Treatment == cd$Reference) %>%
        pull(all_of(metric))
      
      es <- cohen_pooled(x, y)
      
      cohen_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        Family = cd$Family,
        Contrast_order = cd$Contrast_order,
        Contrast = cd$Contrast,
        Target = cd$Target,
        Reference = cd$Reference,
        Context = cd$Context,
        N_target = es["n_target"],
        N_reference = es["n_reference"],
        Mean_target = es["mean_target"],
        Mean_reference = es["mean_reference"],
        Cohen_d_pooled = es["d"],
        Hedges_g = es["g"]
      )
      
      idx <- idx + 1
    }
  }
}

cohen_out <- bind_rows(cohen_rows)

write_out(
  cohen_out,
  "07_Cohen_d_Hedges_g_effect_specific.csv"
)


# Secondary single-treatment diagnostic: descriptive Cohen's d/raw mean
# difference versus the categorical-model W/H/D/C-vs-Ctrl contrast. This is
# retained for transparency but is not the primary Fig. 2a summary in v3.
fig2_effect_diagnostic <- fig2_single4 %>%
  left_join(
    cohen_out %>%
      filter(Contrast %in% c("W - Ctrl", "H - Ctrl", "D - Ctrl", "C - Ctrl")) %>%
      select(
        Group, Metric, Contrast,
        Mean_target, Mean_reference,
        Cohen_d_pooled, Hedges_g
      ),
    by = c("Group", "Metric", "Contrast")
  ) %>%
  mutate(
    Raw_mean_difference = Mean_target - Mean_reference,
    Model_minus_raw_difference = Estimate - Raw_mean_difference,
    Same_direction_raw_and_model =
      sign(Raw_mean_difference) == sign(Estimate),
    Abs_model_to_raw_ratio = ifelse(
      is.finite(Raw_mean_difference) & abs(Raw_mean_difference) > 1e-12,
      abs(Estimate / Raw_mean_difference),
      NA_real_
    )
  )

write_out(
  fig2_effect_diagnostic,
  "07_single_treatment_Cohen_d_vs_LMM_contrast_diagnostic.csv"
)

# ------------------------------------------------------------------------------
# 10B. OVERALL COHEN'S d / HEDGES' g FOR FIG. 2a
#
# The original context-specific Cohen's d values above compare individual
# treatment cells (for example W vs Ctrl). In v3, Fig. 2a instead summarizes the
# OVERALL marginal effect of each global-change factor:
#
#   Warming: Warmed vs Unwarmed, pooling equally represented precipitation x
#            clipping contexts (24 vs 24 observations in the balanced design).
#
#   Half precipitation: Half vs Ambient, pooling equally represented warming x
#                       clipping contexts (16 vs 16 observations).
#
#   Double precipitation: Double vs Ambient, pooling equally represented warming
#                         x clipping contexts (16 vs 16 observations).
#
#   Clipping: Clipped vs Unclipped, pooling equally represented warming x
#             precipitation contexts (24 vs 24 observations).
#
# These pooled-SD Cohen's d / Hedges' g values are DESCRIPTIVE standardized
# magnitudes only. The observations are not treated as independent replicates for
# formal inference. Statistical inference is supplied by the categorical LMM
# overall contrasts already stored in primary_overall.
# ------------------------------------------------------------------------------

overall_effect_def <- tibble(
  Family = c(
    "Warming",
    "Half_precipitation",
    "Double_precipitation",
    "Clipping"
  ),
  Effect = c(
    "Overall_Warming",
    "Overall_Half",
    "Overall_Double",
    "Overall_Clipping"
  ),
  Fig2_treatment = c("W", "H", "D", "C"),
  Target = c("Warmed", "Half", "Double", "Clipped"),
  Reference = c("Unwarmed", "Ambient", "Ambient", "Unclipped"),
  Averaged_over = c(
    "Precipitation x Clipping",
    "Warming x Clipping",
    "Warming x Clipping",
    "Warming x Precipitation"
  )
)

extract_overall_cohen_hedges <- function(d, metric, group_name) {
  dd <- d %>%
    filter(is.finite(.data[[metric]])) %>%
    droplevels()

  get_es_row <- function(family) {
    if (family == "Warming") {
      x <- dd %>% filter(WarmF == "Warmed") %>% pull(all_of(metric))
      y <- dd %>% filter(WarmF == "Unwarmed") %>% pull(all_of(metric))
    } else if (family == "Half_precipitation") {
      x <- dd %>% filter(PrecipF == "Half") %>% pull(all_of(metric))
      y <- dd %>% filter(PrecipF == "Ambient") %>% pull(all_of(metric))
    } else if (family == "Double_precipitation") {
      x <- dd %>% filter(PrecipF == "Double") %>% pull(all_of(metric))
      y <- dd %>% filter(PrecipF == "Ambient") %>% pull(all_of(metric))
    } else if (family == "Clipping") {
      x <- dd %>% filter(ClipF == "Clipped") %>% pull(all_of(metric))
      y <- dd %>% filter(ClipF == "Unclipped") %>% pull(all_of(metric))
    } else {
      stop("Unknown overall effect family: ", family)
    }

    es <- cohen_pooled(x, y)
    def <- overall_effect_def %>% filter(Family == family)

    tibble(
      Group = group_name,
      Metric = metric,
      Family = family,
      Effect = def$Effect[1],
      Fig2_treatment = def$Fig2_treatment[1],
      Target = def$Target[1],
      Reference = def$Reference[1],
      Averaged_over = def$Averaged_over[1],
      N_target = unname(es["n_target"]),
      N_reference = unname(es["n_reference"]),
      Mean_target = unname(es["mean_target"]),
      Mean_reference = unname(es["mean_reference"]),
      SD_target = unname(es["sd_target"]),
      SD_reference = unname(es["sd_reference"]),
      SD_pooled = unname(es["sd_pooled"]),
      Raw_marginal_difference =
        unname(es["mean_target"] - es["mean_reference"]),
      Cohen_d_overall = unname(es["d"]),
      Hedges_g_overall = unname(es["g"])
    )
  }

  bind_rows(lapply(overall_effect_def$Family, get_es_row))
}

overall_cohen_rows <- list()
idx_overall_cohen <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  dg <- dat %>% filter(Group == g) %>% droplevels()

  for (metric in MAIN_METRICS) {
    if (!metric %in% names(dg)) next

    rr <- tryCatch(
      extract_overall_cohen_hedges(dg, metric, g),
      error = function(e) e
    )

    if (!inherits(rr, "error")) {
      overall_cohen_rows[[idx_overall_cohen]] <- rr
      idx_overall_cohen <- idx_overall_cohen + 1
    } else {
      warning(
        "Overall Cohen/Hedges calculation failed for ",
        g, " / ", metric, ": ", conditionMessage(rr)
      )
    }
  }
}

overall_effect_sizes <- bind_rows(overall_cohen_rows)

write_out(
  overall_effect_sizes,
  "07_OVERALL_Cohen_d_Hedges_g.csv"
)

# ------------------------------------------------------------------------------
# 10C. FINAL FIG. 2a OVERALL EFFECT TABLE
#
# primary_overall was generated above from the SAME categorical hierarchical LMM
# used for all primary inference. Each row is the equal-weight average of the
# corresponding context-specific simple effects:
#
# Warming:
#   mean(W-Ctrl, WH-H, WD-D, WC-C, WHC-HC, WDC-DC)
# Half:
#   mean(H-Ctrl, WH-W, HC-C, WHC-WC)
# Double:
#   mean(D-Ctrl, WD-W, DC-C, WDC-WC)
# Clipping:
#   mean(C-Ctrl, WC-W, HC-H, DC-D, WHC-WH, WDC-WD)
#
# P_Holm_four_overall is the Holm correction across these FOUR overall effects
# within each Group x Metric. Therefore, the Fig. 2a bar magnitude and its star
# answer the same overall-factor question, while magnitude and inference remain
# conceptually distinct (descriptive Cohen's d vs hierarchical LMM test).
# ------------------------------------------------------------------------------

fig2a_overall <- primary_overall %>%
  left_join(
    overall_effect_sizes,
    by = c("Group", "Metric", "Family")
  ) %>%
  mutate(
    Effect = case_when(
      Family == "Warming" ~ "Overall_Warming",
      Family == "Half_precipitation" ~ "Overall_Half",
      Family == "Double_precipitation" ~ "Overall_Double",
      Family == "Clipping" ~ "Overall_Clipping",
      TRUE ~ NA_character_
    ),
    Fig2_treatment = case_when(
      Family == "Warming" ~ "W",
      Family == "Half_precipitation" ~ "H",
      Family == "Double_precipitation" ~ "D",
      Family == "Clipping" ~ "C",
      TRUE ~ NA_character_
    ),
    LMM_minus_raw_marginal = Estimate - Raw_marginal_difference,
    Same_direction_raw_and_LMM = case_when(
      !is.finite(Raw_marginal_difference) | !is.finite(Estimate) ~ NA,
      abs(Raw_marginal_difference) < 1e-12 | abs(Estimate) < 1e-12 ~ TRUE,
      TRUE ~ sign(Raw_marginal_difference) == sign(Estimate)
    ),
    Fig2_order = match(Fig2_treatment, c("W", "H", "D", "C"))
  ) %>%
  arrange(Metric, Group, Fig2_order)

write_out(
  fig2a_overall,
  "08_FIG2a_OVERALL_W_H_D_C_effectsize_and_LMM.csv"
)

# Compact Beta_raw/Beta_detrended table for manuscript/figure checking.
fig2a_overall_beta <- fig2a_overall %>%
  filter(Metric %in% c("Beta_raw", "Beta_detrended")) %>%
  select(
    Group, Metric, Family, Effect, Fig2_treatment,
    Target, Reference, Averaged_over,
    N_target, N_reference,
    Mean_target, Mean_reference,
    Raw_marginal_difference,
    SD_target, SD_reference, SD_pooled,
    Cohen_d_overall, Hedges_g_overall,
    Estimate, SE, df, CI_low, CI_high, t_ratio,
    P_raw, P_Holm_four_overall, Significance_Holm_four_overall,
    LMM_minus_raw_marginal, Same_direction_raw_and_LMM
  )

write_out(
  fig2a_overall_beta,
  "08_FIG2a_OVERALL_beta_compact.csv"
)

# ------------------------------------------------------------------------------
# 11. RESPONSE RATIOS
#
# Two forms are output:
# A. Conventional ratio of treatment/reference arithmetic means.
# B. Block-matched log ratio, calculated from the four block-matched pairs.
#
# The block-matched log ratio is closer to the user's original block-paired
# sensitivity workflow.
#
# IMPORTANT:
# RR requires positive values. LogVar_at_ref_raw is therefore back-transformed
# to Var_at_ref_raw before RR calculation.
# ------------------------------------------------------------------------------

rr_source_metric <- function(metric) {
  if (metric == "LogVar_at_ref_raw") {
    "Var_at_ref_raw"
  } else {
    metric
  }
}

unpaired_rr <- function(x, y) {
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  n1 <- length(x)
  n0 <- length(y)
  
  if (
    n1 < 2 || n0 < 2 ||
    any(x <= 0) || any(y <= 0) ||
    mean(x) <= 0 || mean(y) <= 0
  ) {
    return(list(
      available = FALSE,
      reason = "RR requires strictly positive values.",
      RR = NA_real_,
      lnRR = NA_real_,
      SE_lnRR = NA_real_,
      CI_low_RR = NA_real_,
      CI_high_RR = NA_real_,
      P_approx = NA_real_
    ))
  }
  
  mx <- mean(x)
  my <- mean(y)
  
  lnrr <- log(mx / my)
  
  se <- sqrt(
    var(x) / (n1 * mx^2) +
      var(y) / (n0 * my^2)
  )
  
  z <- ifelse(se > 0, lnrr / se, NA_real_)
  p <- ifelse(
    is.finite(z),
    2 * stats::pnorm(abs(z), lower.tail = FALSE),
    NA_real_
  )
  
  list(
    available = TRUE,
    reason = "",
    RR = exp(lnrr),
    lnRR = lnrr,
    SE_lnRR = se,
    CI_low_RR = exp(lnrr - 1.96 * se),
    CI_high_RR = exp(lnrr + 1.96 * se),
    P_approx = p
  )
}

paired_log_rr <- function(d, metric, target, reference) {
  
  pp <- get_matched_pair(
    d,
    metric,
    target,
    reference
  )
  
  if (
    nrow(pp) < 2 ||
    any(pp$Target_value <= 0, na.rm = TRUE) ||
    any(pp$Reference_value <= 0, na.rm = TRUE)
  ) {
    return(list(
      available = FALSE,
      reason = "Paired lnRR requires strictly positive matched values.",
      N = nrow(pp),
      Mean_lnRR = NA_real_,
      SE_lnRR = NA_real_,
      RR_geomean = NA_real_,
      CI_low_RR = NA_real_,
      CI_high_RR = NA_real_,
      P = NA_real_
    ))
  }
  
  lr <- log(
    pp$Target_value /
      pp$Reference_value
  )
  
  lr <- lr[is.finite(lr)]
  
  if (length(lr) < 2 || sd(lr) == 0) {
    return(list(
      available = FALSE,
      reason = "Insufficient variation in paired log ratios.",
      N = length(lr),
      Mean_lnRR = ifelse(length(lr) > 0, mean(lr), NA_real_),
      SE_lnRR = NA_real_,
      RR_geomean = ifelse(length(lr) > 0, exp(mean(lr)), NA_real_),
      CI_low_RR = NA_real_,
      CI_high_RR = NA_real_,
      P = NA_real_
    ))
  }
  
  tt <- t.test(lr, mu = 0)
  
  list(
    available = TRUE,
    reason = "",
    N = length(lr),
    Mean_lnRR = mean(lr),
    SE_lnRR = sd(lr) / sqrt(length(lr)),
    RR_geomean = exp(mean(lr)),
    CI_low_RR = exp(tt$conf.int[1]),
    CI_high_RR = exp(tt$conf.int[2]),
    P = tt$p.value
  )
}

rr_rows <- list()
idx <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  
  dg <- dat %>% filter(Group == g)
  
  for (metric in MAIN_METRICS) {
    
    rr_metric <- rr_source_metric(metric)
    
    if (!rr_metric %in% names(dg)) next
    
    for (j in seq_len(nrow(contrast_def))) {
      
      cd <- contrast_def[j, ]
      
      x <- dg %>%
        filter(Treatment == cd$Target) %>%
        pull(all_of(rr_metric))
      
      y <- dg %>%
        filter(Treatment == cd$Reference) %>%
        pull(all_of(rr_metric))
      
      ur <- unpaired_rr(x, y)
      
      pr <- paired_log_rr(
        dg,
        rr_metric,
        cd$Target,
        cd$Reference
      )
      
      rr_rows[[idx]] <- tibble(
        Group = g,
        Metric = metric,
        RR_metric = rr_metric,
        Family = cd$Family,
        Contrast_order = cd$Contrast_order,
        Contrast = cd$Contrast,
        Target = cd$Target,
        Reference = cd$Reference,
        Context = cd$Context,
        
        RR_available = ur$available,
        RR_reason = ur$reason,
        RR_arithmetic_means = ur$RR,
        lnRR_arithmetic_means = ur$lnRR,
        SE_lnRR_unpaired = ur$SE_lnRR,
        CI_low_RR_unpaired = ur$CI_low_RR,
        CI_high_RR_unpaired = ur$CI_high_RR,
        P_RR_approx = ur$P_approx,
        
        Paired_RR_available = pr$available,
        Paired_RR_reason = pr$reason,
        N_blockmatched = pr$N,
        Mean_blockmatched_lnRR = pr$Mean_lnRR,
        SE_blockmatched_lnRR = pr$SE_lnRR,
        RR_blockmatched_geomean = pr$RR_geomean,
        CI_low_RR_blockmatched = pr$CI_low_RR,
        CI_high_RR_blockmatched = pr$CI_high_RR,
        P_RR_blockmatched = pr$P
      )
      
      idx <- idx + 1
    }
  }
}

rr_out <- bind_rows(rr_rows) %>%
  group_by(Group, Metric, Family) %>%
  mutate(
    P_RR_blockmatched_Holm =
      p.adjust(P_RR_blockmatched, method = "holm")
  ) %>%
  ungroup()

write_out(
  rr_out,
  "08_response_ratio_effect_specific.csv"
)

# ------------------------------------------------------------------------------
# 12. OPTIONAL LEGACY ieggr::compv
#
# Reproduces the style of the user's original code but now uses the correct
# effect-specific comparison pairs instead of comparing everything with Ctrl.
# ------------------------------------------------------------------------------

legacy_rows <- list()

if (
  RUN_LEGACY_IEGGR &&
  requireNamespace("ieggr", quietly = TRUE)
) {
  
  compv_fun <- tryCatch(
    get("compv", envir = asNamespace("ieggr")),
    error = function(e) NULL
  )
  
  if (!is.null(compv_fun)) {
    
    idx <- 1
    
    for (g in unique(as.character(na.omit(dat$Group)))) {
      
      dg <- dat %>% filter(Group == g)
      
      # Use only metrics for which direct RR is conceptually meaningful.
      legacy_metrics <- intersect(
        c(
          "Beta_raw",
          "Beta_detrended",
          "CTV_Hellinger_RMS",
          "WeightedTaxonCV_raw",
          "MedianTaxonCV_raw",
          "MeanPairwise_Hellinger",
          "MeanAnnual_Hellinger",
          "MeanPairwise_BrayCurtis",
          "MeanAnnual_BrayCurtis"
        ),
        names(dg)
      )
      
      for (metric in legacy_metrics) {
        
        for (j in seq_len(nrow(contrast_def))) {
          
          cd <- contrast_def[j, ]
          
          tx <- dg %>%
            filter(Treatment == cd$Target) %>%
            arrange(block)
          
          rf <- dg %>%
            filter(Treatment == cd$Reference) %>%
            arrange(block)
          
          # Unpaired legacy compv.
          cu <- tryCatch(
            as.data.frame(
              compv_fun(
                tx[[metric]],
                rf[[metric]],
                paired = FALSE,
                alternative = c(
                  "two.sided", "less", "greater"
                ),
                method = c("t", "w", "d", "r"),
                nworker = 1,
                rr.fill.zero = c("yes", "no")
              )
            ),
            error = function(e) {
              data.frame(
                legacy_error = conditionMessage(e)
              )
            }
          )
          
          cu$Group <- g
          cu$Metric <- metric
          cu$Family <- cd$Family
          cu$Contrast <- cd$Contrast
          cu$Legacy_mode <- "unpaired"
          
          legacy_rows[[idx]] <- cu
          idx <- idx + 1
          
          # Block-matched legacy compv.
          cp <- tryCatch(
            as.data.frame(
              compv_fun(
                tx[[metric]],
                rf[[metric]],
                paired = TRUE,
                alternative = c(
                  "two.sided", "less", "greater"
                ),
                method = c("t", "w", "d", "r"),
                nworker = 1,
                rr.fill.zero = c("yes", "no")
              )
            ),
            error = function(e) {
              data.frame(
                legacy_error = conditionMessage(e)
              )
            }
          )
          
          cp$Group <- g
          cp$Metric <- metric
          cp$Family <- cd$Family
          cp$Contrast <- cd$Contrast
          cp$Legacy_mode <- "block_matched"
          
          legacy_rows[[idx]] <- cp
          idx <- idx + 1
        }
      }
    }
  }
}

if (length(legacy_rows) > 0) {
  
  # bind_rows tolerates version-dependent compv column sets.
  legacy_out <- bind_rows(legacy_rows)
  
  write_out(
    legacy_out,
    "09_legacy_ieggr_compv_effect_specific.csv"
  )
  
} else {
  
  writeLines(
    c(
      "Legacy ieggr::compv analysis was not run.",
      "Either RUN_LEGACY_IEGGR was FALSE or package 'ieggr' was unavailable.",
      "This does not affect the primary hierarchical analysis."
    ),
    file.path(
      OUTDIR,
      "09_legacy_ieggr_compv_NOT_RUN.txt"
    )
  )
}

# ------------------------------------------------------------------------------
# 13. BETA-SE PROPAGATION WITH MULTILEVEL META-REGRESSION
# ------------------------------------------------------------------------------
# The uncertainty of each subplot-level beta estimate is propagated through
# V = SE_Beta^2 in rma.mv(). The fixed-effect structure mirrors the categorical
# primary model so that Half and Double precipitation are estimated independently.
# This remains a SENSITIVITY ANALYSIS and does not replace the primary LMM.
# ------------------------------------------------------------------------------

BETA_SE_PAIRS <- BETA_SE_PAIRS[
  names(BETA_SE_PAIRS) %in% names(dat) &
    unname(BETA_SE_PAIRS) %in% names(dat)
]

if (length(BETA_SE_PAIRS) == 0) {
  warning(
    "No usable Beta + SE_Beta pair was found in the input table. ",
    "The Beta-SE propagation sensitivity section will be skipped."
  )
}

fit_beta_se_model <- function(d, beta_col, se_col) {
  dd <- d %>%
    filter(
      is.finite(.data[[beta_col]]),
      is.finite(.data[[se_col]]),
      .data[[se_col]] > 0
    ) %>%
    droplevels()
  
  if (nrow(dd) < 20) {
    return(list(fit = NULL, data = dd, error = "Too few observations."))
  }
  
  fit <- tryCatch(
    metafor::rma.mv(
      yi = dd[[beta_col]],
      V = dd[[se_col]]^2,
      mods = ~ WarmF * PrecipF * ClipF,
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
    return(list(fit = NULL, data = dd, error = conditionMessage(fit)))
  }
  list(fit = fit, data = dd, error = NA_character_)
}

# One-row factor representation of a manuscript treatment label. Factor levels
# are fixed explicitly so model.matrix() reproduces the primary design matrix.
treatment_factor_row <- function(label) {
  grid <- tibble(
    Treatment = TREATMENT_LEVELS,
    WarmF = factor(
      c("Unwarmed", "Warmed", "Unwarmed", "Unwarmed", "Unwarmed",
        "Warmed", "Warmed", "Warmed", "Unwarmed", "Unwarmed", "Warmed", "Warmed"),
      levels = c("Unwarmed", "Warmed")
    ),
    PrecipF = factor(
      c("Ambient", "Ambient", "Half", "Double", "Ambient",
        "Half", "Double", "Ambient", "Half", "Double", "Half", "Double"),
      levels = c("Ambient", "Half", "Double")
    ),
    ClipF = factor(
      c("Unclipped", "Unclipped", "Unclipped", "Unclipped", "Clipped",
        "Unclipped", "Unclipped", "Clipped", "Clipped", "Clipped", "Clipped", "Clipped"),
      levels = c("Unclipped", "Clipped")
    )
  )
  out <- grid %>% filter(Treatment == label)
  if (nrow(out) != 1) stop("Unknown treatment label: ", label)
  out
}

meta_pair_weight <- function(fit, target, reference) {
  tt <- treatment_factor_row(target)
  rr <- treatment_factor_row(reference)
  mm <- model.matrix(
    ~ WarmF * PrecipF * ClipF,
    data = bind_rows(tt, rr)
  )
  drow <- mm[1, , drop = TRUE] - mm[2, , drop = TRUE]
  
  bn <- names(coef(fit))
  L <- rep(0, length(bn))
  names(L) <- bn
  
  for (nm in bn) {
    source_name <- if (nm == "intrcpt") "(Intercept)" else nm
    if (source_name %in% names(drow)) {
      L[nm] <- drow[source_name]
    } else {
      stop(
        "Could not map meta-regression coefficient '", nm,
        "' to the categorical WarmF*PrecipF*ClipF design matrix."
      )
    }
  }
  L
}

meta_model_contrast <- function(fit, L, label) {
  Lmat <- matrix(
    as.numeric(L),
    nrow = 1,
    dimnames = list(label, names(L))
  )
  aa <- tryCatch(stats::anova(fit, X = Lmat), error = function(e) e)
  bb <- tryCatch(as.numeric(stats::coef(fit)), error = function(e) NULL)
  vv <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
  
  est <- se <- NA_real_
  if (!is.null(bb) && length(bb) == length(L)) {
    est <- sum(as.numeric(L) * bb)
  }
  if (!is.null(vv) && nrow(vv) == length(L) && ncol(vv) == length(L)) {
    se2 <- as.numeric(t(as.numeric(L)) %*% vv %*% as.numeric(L))
    if (is.finite(se2) && se2 >= 0) se <- sqrt(se2)
  }
  
  p <- ddf <- NA_real_
  test_error <- NA_character_
  if (inherits(aa, "error")) {
    test_error <- conditionMessage(aa)
  } else {
    if (!is.null(aa$pval)) p <- as.numeric(aa$pval[1])
    if (!is.null(aa$ddf)) ddf <- as.numeric(aa$ddf[1])
    else if (!is.null(aa$dfs)) ddf <- as.numeric(aa$dfs[1])
    else if (!is.null(aa$QMdf) && length(aa$QMdf) >= 2) ddf <- as.numeric(aa$QMdf[2])
  }
  
  if (!is.finite(p) && is.finite(est) && is.finite(se) && se > 0) {
    stat <- est / se
    p <- if (is.finite(ddf) && ddf > 0) {
      2 * pt(abs(stat), df = ddf, lower.tail = FALSE)
    } else {
      2 * pnorm(abs(stat), lower.tail = FALSE)
    }
  }
  
  crit <- if (is.finite(ddf) && ddf > 0) qt(0.975, df = ddf) else qnorm(0.975)
  lo <- if (is.finite(est) && is.finite(se)) est - crit * se else NA_real_
  hi <- if (is.finite(est) && is.finite(se)) est + crit * se else NA_real_
  
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
idx_diag <- idx_con <- idx_over <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  dg <- dat %>% filter(Group == g)
  
  for (beta_metric in names(BETA_SE_PAIRS)) {
    se_metric <- unname(BETA_SE_PAIRS[[beta_metric]])
    ans <- fit_beta_se_model(dg, beta_metric, se_metric)
    
    if (is.null(ans$fit)) {
      meta_diag_rows[[idx_diag]] <- tibble(
        Group = g,
        Metric = beta_metric,
        SE_metric = se_metric,
        N = nrow(ans$data),
        N_blocks = n_distinct(ans$data$Block),
        N_wholeplots = n_distinct(ans$data$BlockWholePlotID),
        Median_input_SE = if (nrow(ans$data) > 0) median(ans$data[[se_metric]], na.rm = TRUE) else NA_real_,
        Mean_input_SE = if (nrow(ans$data) > 0) mean(ans$data[[se_metric]], na.rm = TRUE) else NA_real_,
        Empirical_SD_Beta = if (nrow(ans$data) > 1) sd(ans$data[[beta_metric]], na.rm = TRUE) else NA_real_,
        Block_variance = NA_real_,
        WholePlot_within_Block_variance = NA_real_,
        Plot_heterogeneity_variance = NA_real_,
        QE = NA_real_,
        QE_P = NA_real_,
        Inference_method = "rma.mv; V=SE_Beta^2; categorical treatment factors; test=t; dfs=contain",
        Error = ans$error
      )
      idx_diag <- idx_diag + 1
      next
    }
    
    fit <- ans$fit
    dd <- ans$data
    sig2 <- fit$sigma2
    meta_diag_rows[[idx_diag]] <- tibble(
      Group = g,
      Metric = beta_metric,
      SE_metric = se_metric,
      N = nrow(dd),
      N_blocks = n_distinct(dd$Block),
      N_wholeplots = n_distinct(dd$BlockWholePlotID),
      Median_input_SE = median(dd[[se_metric]], na.rm = TRUE),
      Mean_input_SE = mean(dd[[se_metric]], na.rm = TRUE),
      Empirical_SD_Beta = sd(dd[[beta_metric]], na.rm = TRUE),
      Block_variance = ifelse(length(sig2) >= 1, sig2[1], NA),
      WholePlot_within_Block_variance = ifelse(length(sig2) >= 2, sig2[2], NA),
      Plot_heterogeneity_variance = ifelse(length(sig2) >= 3, sig2[3], NA),
      QE = fit$QE,
      QE_P = fit$QEp,
      Inference_method = "rma.mv; V=SE_Beta^2; categorical treatment factors; test=t; dfs=contain",
      Error = NA_character_
    )
    idx_diag <- idx_diag + 1
    
    family_weights <- list()
    for (fam in unique(contrast_def$Family)) {
      defs <- contrast_def %>% filter(Family == fam) %>% arrange(Contrast_order)
      fw <- list()
      for (j in seq_len(nrow(defs))) {
        L <- meta_pair_weight(fit, defs$Target[j], defs$Reference[j])
        fw[[j]] <- L
        rr <- meta_model_contrast(fit, L, defs$Contrast[j]) %>%
          mutate(
            Group = g,
            Metric = beta_metric,
            Family = fam,
            Contrast_order = defs$Contrast_order[j],
            Target = defs$Target[j],
            Reference = defs$Reference[j],
            Context = defs$Context[j],
            .before = 1
          )
        meta_contrast_rows[[idx_con]] <- rr
        idx_con <- idx_con + 1
      }
      family_weights[[fam]] <- fw
      Lavg <- make_average_weight(fw)
      oo <- meta_model_contrast(fit, Lavg, paste0("Overall_", fam)) %>%
        mutate(Group = g, Metric = beta_metric, Family = fam, .before = 1)
      meta_overall_rows[[idx_over]] <- oo
      idx_over <- idx_over + 1
    }
  }
}

meta_diag <- bind_rows(meta_diag_rows)
meta_contrasts <- bind_rows(meta_contrast_rows)
meta_overall <- bind_rows(meta_overall_rows)

if (nrow(meta_contrasts) > 0) {
  meta_contrasts <- meta_contrasts %>%
    group_by(Group, Metric, Family) %>%
    mutate(
      P_SEprop_Holm_within_family = p.adjust(P_SEprop, method = "holm"),
      Significance_SEprop_Holm = sig_label(P_SEprop_Holm_within_family)
    ) %>%
    ungroup() %>%
    group_by(Group, Metric) %>%
    mutate(P_SEprop_Holm_all20 = p.adjust(P_SEprop, method = "holm")) %>%
    ungroup()
}

if (nrow(meta_overall) > 0) {
  meta_overall <- meta_overall %>%
    group_by(Group, Metric) %>%
    mutate(
      P_SEprop_Holm_four_overall = p.adjust(P_SEprop, method = "holm"),
      Significance_SEprop_Holm = sig_label(P_SEprop_Holm_four_overall)
    ) %>%
    ungroup()
}

write_out(meta_diag, "10_Beta_SE_propagation_model_diagnostics.csv")
write_out(meta_contrasts, "10_Beta_SE_propagation_effect_specific_contrasts.csv")
write_out(meta_overall, "10_Beta_SE_propagation_overall_directional_effects.csv")

old_cr2_files <- file.path(
  OUTDIR,
  c(
    "10_Beta_SE_CR2_model_diagnostics.csv",
    "10_Beta_SE_CR2_effect_specific_contrasts.csv",
    "10_Beta_SE_CR2_overall_directional_effects.csv"
  )
)
unlink(old_cr2_files[file.exists(old_cr2_files)])

writeLines(
  c(
    "CR2 was intentionally not used in the final Beta-SE propagation analysis.",
    "Beta-estimation uncertainty is propagated directly through V = SE_Beta^2 in rma.mv().",
    "The fixed-effect structure mirrors the categorical primary model (WarmF * PrecipF * ClipF).",
    "With Block as the highest-level random effect there are only four independent top-level clusters,",
    "so CR2 cluster-robust inference is not appropriate for the manuscript's standard nested model."
  ),
  file.path(OUTDIR, "10_Beta_SE_CR2_NOT_USED.txt")
)

# ------------------------------------------------------------------------------
# 14. UNIFIED EFFECT-SPECIFIC TABLE
# ------------------------------------------------------------------------------

unified <- primary_contrasts %>%
  left_join(
    cohen_out %>%
      select(
        Group, Metric, Family, Contrast,
        Cohen_d_pooled, Hedges_g,
        Mean_target, Mean_reference
      ),
    by = c(
      "Group", "Metric", "Family", "Contrast"
    )
  ) %>%
  left_join(
    blockmatched %>%
      select(
        Group, Metric, Family, Contrast,
        Mean_block_difference,
        SE_block_difference,
        P_blockmatched,
        P_blockmatched_Holm,
        Cohen_dz_blockmatched,
        All_same_wholeplot
      ),
    by = c(
      "Group", "Metric", "Family", "Contrast"
    )
  ) %>%
  left_join(
    rr_out %>%
      select(
        Group, Metric, Family, Contrast,
        RR_metric,
        RR_arithmetic_means,
        P_RR_approx,
        RR_blockmatched_geomean,
        P_RR_blockmatched,
        P_RR_blockmatched_Holm
      ),
    by = c(
      "Group", "Metric", "Family", "Contrast"
    )
  )

if (nrow(meta_contrasts) > 0) {
  
  unified <- unified %>%
    left_join(
      meta_contrasts %>%
        select(
          Group, Metric, Family, Contrast,
          Estimate_SEprop,
          SE_SEprop,
          CI_low_SEprop,
          CI_high_SEprop,
          df_SEprop,
          P_SEprop,
          P_SEprop_Holm_within_family
        ),
      by = c(
        "Group", "Metric", "Family", "Contrast"
      )
    )
}

write_out(
  unified,
  "11_UNIFIED_effect_specific_all_evidence.csv"
)

# Unified overall effects.
# Add descriptive overall Cohen's d / Hedges' g so the reviewer/manuscript table
# places standardized magnitude and hierarchical inference side by side.
unified_overall <- primary_overall %>%
  left_join(
    overall_effect_sizes %>%
      select(
        Group, Metric, Family, Effect, Fig2_treatment,
        Target, Reference, Averaged_over,
        N_target, N_reference, Mean_target, Mean_reference,
        Raw_marginal_difference, SD_pooled,
        Cohen_d_overall, Hedges_g_overall
      ),
    by = c("Group", "Metric", "Family")
  ) %>%
  left_join(
    block_overall,
    by = c(
      "Group", "Metric", "Family"
    ),
    suffix = c(
      "_hierarchical",
      "_blockmatched"
    )
  )

if (nrow(meta_overall) > 0) {
  
  unified_overall <- unified_overall %>%
    left_join(
      meta_overall %>%
        select(
          Group, Metric, Family,
          Estimate_SEprop,
          SE_SEprop,
          CI_low_SEprop,
          CI_high_SEprop,
          df_SEprop,
          P_SEprop,
          P_SEprop_Holm_four_overall
        ),
      by = c(
        "Group", "Metric", "Family"
      )
    )
}

write_out(
  unified_overall,
  "12_UNIFIED_overall_directional_effects.csv"
)

# ------------------------------------------------------------------------------
# 15. OPTIONAL BETA-COMMUNITY VARIABILITY ASSOCIATION
#
# This is descriptive/supplementary and should not be used as proof that beta
# equals community variability.
# ------------------------------------------------------------------------------

assoc_rows <- list()
idx <- 1

if (
  all(
    c(
      "Beta_raw",
      "CTV_Hellinger_RMS"
    ) %in% names(dat)
  )
) {
  
  for (g in unique(as.character(na.omit(dat$Group)))) {
    
    dg <- dat %>%
      filter(Group == g) %>%
      mutate(
        Beta_raw_z =
          as.numeric(scale(Beta_raw))
      )
    
    m <- tryCatch(
      lmerTest::lmer(
        CTV_Hellinger_RMS ~
          Beta_raw_z +
          WarmF * PrecipF * ClipF +
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
      
      if ("Beta_raw_z" %in% rownames(co)) {
        
        rr <- co["Beta_raw_z", , drop = FALSE]
        
        assoc_rows[[idx]] <- tibble(
          Group = g,
          Estimate_beta_z = rr[1, "Estimate"],
          SE = rr[1, "Std. Error"],
          df = rr[1, "df"],
          t = rr[1, "t value"],
          P = rr[1, "Pr(>|t|)"],
          Singular =
            lme4::isSingular(m, tol = 1e-6)
        )
        
        idx <- idx + 1
      }
    }
  }
}

assoc_out <- bind_rows(assoc_rows)

if (nrow(assoc_out) > 0) {
  write_out(
    assoc_out,
    "13_beta_CTV_association_adjusted_descriptive.csv"
  )
}

# ------------------------------------------------------------------------------
# 15B. BETA / REFERENCE-VARIANCE CO-VARIATION UNDER WARMING
#
# Reviewer question addressed here:
#   "Does the intercept and slope of TPL co-vary under warming?"
#
# IMPORTANT INTERPRETATION
# ------------------------
# This section does NOT use the conventional raw TPL intercept, because that
# intercept corresponds to the fitted variance at mean relative abundance = 1,
# far outside the abundance range of most taxa and is strongly parameterization
# dependent. Instead, it uses LogVar_at_ref_raw, the fitted log10 temporal
# variance at the group-specific reference abundance already used elsewhere in
# the revision.
#
# The primary co-variation question is phrased as a change-change association:
# within each microbial group, block, precipitation level and clipping level,
# calculate the warming-induced change in Beta_raw and the corresponding
# warming-induced change in LogVar_at_ref_raw. Then test whether larger
# warming-induced changes in Beta_raw are associated with larger changes in the
# reference-abundance variance.
#
# There are 4 blocks x 3 precipitation levels x 2 clipping levels = 24 paired
# warming contrasts per microbial group. The unclipped and clipped contrasts at
# the same Block x Precipitation level originate from the same pair of physical
# whole plots, so PairID is included as a random intercept. Block is treated as
# fixed in this reviewer-specific supplementary model because there are only four
# blocks and the target of inference here is the Delta_Beta association rather
# than a treatment main effect.
#
# Primary reviewer-specific model:
#   Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID)
#
# Null hypothesis for co-variation:
#   H0: coefficient(Delta_Beta) = 0
#
# A block-fixed OLS model without the PairID random intercept is also reported as
# a simple sensitivity check. Pearson and Spearman correlations are descriptive
# only; the adjusted mixed model is the preferred inferential result.
# ------------------------------------------------------------------------------

covary_required <- c(
  "Beta_raw",
  "LogVar_at_ref_raw"
)

if (all(covary_required %in% names(dat))) {
  # ---------------------------------------------------------------------------
  # 15B.1 Pairing QC before constructing warming changes
  # ---------------------------------------------------------------------------
  covary_cell_qc <- dat %>%
    filter(
      is.finite(Beta_raw),
      is.finite(LogVar_at_ref_raw)
    ) %>%
    count(
      Group, Block,
      Precipitation, Clipping, Warming,
      name = "N"
    )
  
  write_out(
    covary_cell_qc,
    "13b_Beta_LogVar_covary_cell_QC.csv"
  )
  
  covary_bad_cells <- covary_cell_qc %>%
    filter(N != 1)
  
  if (nrow(covary_bad_cells) > 0) {
    warning(
      "Beta/LogVar co-variation pairing found non-unique cells. ",
      "Inspect 13b_Beta_LogVar_covary_cell_QC.csv."
    )
  }
  
  # ---------------------------------------------------------------------------
  # 15B.2 Construct block-matched warming-induced changes
  # ---------------------------------------------------------------------------
  covary_unwarmed <- dat %>%
    filter(
      Warming == 0,
      is.finite(Beta_raw),
      is.finite(LogVar_at_ref_raw)
    ) %>%
    transmute(
      Group,
      Block,
      Precipitation,
      Clipping,
      PrecipF,
      ClipF,
      Beta_unwarmed = Beta_raw,
      LogVar_ref_unwarmed = LogVar_at_ref_raw,
      WholePlot_unwarmed = as.character(WholePlotID),
      Plot_unwarmed = as.character(Plot)
    )
  
  covary_warmed <- dat %>%
    filter(
      Warming == 1,
      is.finite(Beta_raw),
      is.finite(LogVar_at_ref_raw)
    ) %>%
    transmute(
      Group,
      Block,
      Precipitation,
      Clipping,
      PrecipF,
      ClipF,
      Beta_warmed = Beta_raw,
      LogVar_ref_warmed = LogVar_at_ref_raw,
      WholePlot_warmed = as.character(WholePlotID),
      Plot_warmed = as.character(Plot)
    )
  
  warming_covary_data <- inner_join(
    covary_warmed,
    covary_unwarmed,
    by = c(
      "Group", "Block",
      "Precipitation", "Clipping",
      "PrecipF", "ClipF"
    )
  ) %>%
    mutate(
      Delta_Beta = Beta_warmed - Beta_unwarmed,
      Delta_LogVar_ref =
        LogVar_ref_warmed - LogVar_ref_unwarmed,
      PairID = interaction(
        Block,
        PrecipF,
        drop = TRUE,
        sep = ":P"
      ),
      Precipitation_label = case_when(
        Precipitation == -0.5 ~ "Half",
        Precipitation == 0    ~ "Ambient",
        Precipitation == 1    ~ "Double",
        TRUE ~ NA_character_
      ),
      Clipping_label = ifelse(
        Clipping == 1,
        "Clipped",
        "Unclipped"
      ),
      Context = paste(
        Precipitation_label,
        Clipping_label,
        sep = " | "
      )
    ) %>%
    arrange(Group, Block, Precipitation, Clipping)
  
  write_out(
    warming_covary_data,
    "13b_Beta_LogVar_covary_warming_change_data.csv"
  )
  
  covary_pair_qc <- warming_covary_data %>%
    group_by(Group) %>%
    summarise(
      N_warming_contrasts = n(),
      N_blocks = n_distinct(Block),
      N_pairIDs = n_distinct(PairID),
      N_precipitation_levels =
        n_distinct(Precipitation),
      N_clipping_levels =
        n_distinct(Clipping),
      Min_Delta_Beta = min(Delta_Beta, na.rm = TRUE),
      Max_Delta_Beta = max(Delta_Beta, na.rm = TRUE),
      SD_Delta_Beta = sd(Delta_Beta, na.rm = TRUE),
      Min_Delta_LogVar_ref =
        min(Delta_LogVar_ref, na.rm = TRUE),
      Max_Delta_LogVar_ref =
        max(Delta_LogVar_ref, na.rm = TRUE),
      SD_Delta_LogVar_ref =
        sd(Delta_LogVar_ref, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Expected_24_contrasts =
        N_warming_contrasts == 24,
      Expected_12_pairIDs =
        N_pairIDs == 12
    )
  
  write_out(
    covary_pair_qc,
    "13b_Beta_LogVar_covary_pairing_QC.csv"
  )
  
  # ---------------------------------------------------------------------------
  # 15B.3 Group-specific adjusted co-variation model
  # ---------------------------------------------------------------------------
  fit_covary_group <- function(dg, group_name) {
    dd <- dg %>%
      filter(
        is.finite(Delta_Beta),
        is.finite(Delta_LogVar_ref),
        !is.na(PrecipF),
        !is.na(ClipF)
      ) %>%
      mutate(
        Block = droplevels(factor(Block)),
        PrecipF = droplevels(factor(PrecipF, levels = c("Ambient", "Half", "Double"))),
        ClipF = droplevels(factor(ClipF, levels = c("Unclipped", "Clipped"))),
        PairID = droplevels(factor(PairID))
      )
    
    # Guard against a degenerate predictor.
    if (
      nrow(dd) < 12 ||
      !is.finite(sd(dd$Delta_Beta)) ||
      sd(dd$Delta_Beta) == 0
    ) {
      return(list(
        result = tibble(
          Group = group_name,
          N = nrow(dd),
          N_blocks = n_distinct(dd$Block),
          N_pairIDs = n_distinct(dd$PairID),
          Estimate_DeltaBeta = NA_real_,
          SE = NA_real_,
          df = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_,
          t_ratio = NA_real_,
          P_adjusted = NA_real_,
          PairID_variance = NA_real_,
          Residual_variance = NA_real_,
          Singular = NA,
          Pearson_r = NA_real_,
          Pearson_P = NA_real_,
          Spearman_rho = NA_real_,
          Spearman_P = NA_real_,
          OLS_Estimate_DeltaBeta = NA_real_,
          OLS_SE = NA_real_,
          OLS_df = NA_real_,
          OLS_t = NA_real_,
          OLS_P = NA_real_,
          OLS_R2 = NA_real_,
          OLS_Adj_R2 = NA_real_,
          Error = "Too few observations or zero SD in Delta_Beta."
        ),
        model = NULL,
        ols = NULL,
        data = dd
      ))
    }
    
    form_mixed <-
      Delta_LogVar_ref ~
      Delta_Beta +
      PrecipF * ClipF +
      Block +
      (1 | PairID)
    
    m <- tryCatch(
      lmerTest::lmer(
        form_mixed,
        data = dd,
        REML = TRUE
      ),
      error = function(e) e
    )
    
    # Block-fixed OLS sensitivity without PairID random intercept.
    m_ols <- tryCatch(
      lm(
        Delta_LogVar_ref ~
          Delta_Beta +
          PrecipF * ClipF +
          Block,
        data = dd
      ),
      error = function(e) e
    )
    
    # Descriptive bivariate associations.
    pearson <- tryCatch(
      cor.test(
        dd$Delta_Beta,
        dd$Delta_LogVar_ref,
        method = "pearson"
      ),
      error = function(e) NULL
    )
    
    spearman <- tryCatch(
      cor.test(
        dd$Delta_Beta,
        dd$Delta_LogVar_ref,
        method = "spearman",
        exact = FALSE
      ),
      error = function(e) NULL
    )
    
    if (inherits(m, "error")) {
      mixed_est <- mixed_se <- mixed_df <- mixed_t <- mixed_p <- NA_real_
      mixed_lo <- mixed_hi <- NA_real_
      pair_var <- resid_var <- NA_real_
      singular <- NA
      mixed_error <- conditionMessage(m)
    } else {
      sm <- as.data.frame(summary(m)$coefficients)
      
      if ("Delta_Beta" %in% rownames(sm)) {
        rr <- sm["Delta_Beta", , drop = FALSE]
        mixed_est <- rr[1, "Estimate"]
        mixed_se  <- rr[1, "Std. Error"]
        mixed_df  <- rr[1, "df"]
        mixed_t   <- rr[1, "t value"]
        mixed_p   <- rr[1, "Pr(>|t|)"]
        crit <- qt(0.975, df = mixed_df)
        mixed_lo <- mixed_est - crit * mixed_se
        mixed_hi <- mixed_est + crit * mixed_se
      } else {
        mixed_est <- mixed_se <- mixed_df <- mixed_t <- mixed_p <- NA_real_
        mixed_lo <- mixed_hi <- NA_real_
      }
      
      vc <- as.data.frame(VarCorr(m))
      pair_var <- vc %>%
        filter(grp == "PairID") %>%
        pull(vcov)
      if (length(pair_var) == 0) pair_var <- NA_real_
      
      resid_var <- vc %>%
        filter(grp == "Residual") %>%
        pull(vcov)
      if (length(resid_var) == 0) resid_var <- NA_real_
      
      singular <- lme4::isSingular(m, tol = 1e-6)
      mixed_error <- NA_character_
    }
    
    if (inherits(m_ols, "error")) {
      ols_est <- ols_se <- ols_df <- ols_t <- ols_p <- NA_real_
      ols_r2 <- ols_adj_r2 <- NA_real_
      ols_error <- conditionMessage(m_ols)
    } else {
      so <- summary(m_ols)
      co <- as.data.frame(so$coefficients)
      
      if ("Delta_Beta" %in% rownames(co)) {
        ro <- co["Delta_Beta", , drop = FALSE]
        ols_est <- ro[1, "Estimate"]
        ols_se  <- ro[1, "Std. Error"]
        ols_t   <- ro[1, "t value"]
        ols_p   <- ro[1, "Pr(>|t|)"]
        ols_df  <- df.residual(m_ols)
      } else {
        ols_est <- ols_se <- ols_t <- ols_p <- ols_df <- NA_real_
      }
      
      ols_r2 <- so$r.squared
      ols_adj_r2 <- so$adj.r.squared
      ols_error <- NA_character_
    }
    
    result <- tibble(
      Group = group_name,
      N = nrow(dd),
      N_blocks = n_distinct(dd$Block),
      N_pairIDs = n_distinct(dd$PairID),
      Model =
        "Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID)",
      Null_hypothesis =
        "H0: coefficient(Delta_Beta) = 0",
      Estimate_DeltaBeta = mixed_est,
      SE = mixed_se,
      df = mixed_df,
      CI_low = mixed_lo,
      CI_high = mixed_hi,
      t_ratio = mixed_t,
      P_adjusted = mixed_p,
      PairID_variance = pair_var[1],
      Residual_variance = resid_var[1],
      Singular = singular,
      Pearson_r = if (!is.null(pearson)) unname(pearson$estimate) else NA_real_,
      Pearson_P = if (!is.null(pearson)) pearson$p.value else NA_real_,
      Spearman_rho = if (!is.null(spearman)) unname(spearman$estimate) else NA_real_,
      Spearman_P = if (!is.null(spearman)) spearman$p.value else NA_real_,
      OLS_Estimate_DeltaBeta = ols_est,
      OLS_SE = ols_se,
      OLS_df = ols_df,
      OLS_t = ols_t,
      OLS_P = ols_p,
      OLS_R2 = ols_r2,
      OLS_Adj_R2 = ols_adj_r2,
      Error = paste(
        na.omit(c(
          mixed_error,
          if (!is.na(ols_error)) paste0("OLS: ", ols_error) else NA_character_
        )),
        collapse = " | "
      )
    )
    
    if (identical(result$Error, "")) {
      result$Error <- NA_character_
    }
    
    list(
      result = result,
      model = if (!inherits(m, "error")) m else NULL,
      ols = if (!inherits(m_ols, "error")) m_ols else NULL,
      data = dd
    )
  }
  
  covary_results_list <- list()
  covary_models <- list()
  
  idx_cov <- 1
  
  for (g in unique(as.character(na.omit(warming_covary_data$Group)))) {
    dg <- warming_covary_data %>%
      filter(Group == g)
    
    ans_cov <- fit_covary_group(dg, g)
    
    covary_results_list[[idx_cov]] <- ans_cov$result
    covary_models[[g]] <- ans_cov$model
    
    idx_cov <- idx_cov + 1
  }
  
  covary_results <- bind_rows(covary_results_list) %>%
    mutate(
      P_Holm_across_groups =
        p.adjust(P_adjusted, method = "holm"),
      Significance_raw = sig_label(P_adjusted),
      Significance_Holm =
        sig_label(P_Holm_across_groups),
      Direction = case_when(
        is.na(Estimate_DeltaBeta) ~ NA_character_,
        Estimate_DeltaBeta > 0 ~ "Positive",
        Estimate_DeltaBeta < 0 ~ "Negative",
        TRUE ~ "Zero"
      )
    )
  
  write_out(
    covary_results,
    "13b_Beta_LogVar_covary_adjusted_tests.csv"
  )
  
  # ---------------------------------------------------------------------------
  # 15B.4 Supplementary visualization
  # ---------------------------------------------------------------------------
  p_covary <- ggplot(
    warming_covary_data,
    aes(
      x = Delta_Beta,
      y = Delta_LogVar_ref,
      shape = factor(Clipping)
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted",
      linewidth = 0.35
    ) +
    geom_point(
      size = 2.1,
      alpha = 0.80
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.7
    ) +
    facet_wrap(
      ~ Group,
      scales = "free",
      ncol = 3
    ) +
    scale_shape_manual(
      values = c("0" = 16, "1" = 17),
      labels = c(
        "0" = "Unclipped",
        "1" = "Clipped"
      )
    ) +
    labs(
      x = expression(Delta*beta["warming"]),
      y = expression(Delta*log[10]*Var["ref,warming"]),
      shape = "Clipping"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
  
  ggsave(
    file.path(
      FIGDIR,
      "Fig_S13d_Beta_LogVar_covary_warming_changes.pdf"
    ),
    plot = p_covary,
    width = 8.2,
    height = 5.8
  )
  
  ggsave(
    file.path(
      FIGDIR,
      "Fig_S13d_Beta_LogVar_covary_warming_changes.png"
    ),
    plot = p_covary,
    width = 8.2,
    height = 5.8,
    dpi = 300
  )
  
} else {
  warning(
    "Beta/LogVar co-variation analysis skipped because Beta_raw and/or ",
    "LogVar_at_ref_raw is missing from the input table."
  )
  
  warming_covary_data <- tibble()
  covary_results <- tibble()
}

# ------------------------------------------------------------------------------
# 15C. FIG. 2a OVERALL FACTOR EFFECT DATA/PLOTS
# ------------------------------------------------------------------------------
# PRIMARY Fig. 2a representation in v3:
#   bar height = descriptive overall Cohen's d
#   asterisk   = Holm-adjusted P from the corresponding categorical LMM overall
#                contrast (four overall effects adjusted within Group x Metric)
#
# A second optional plot shows the LMM overall contrast (Delta beta) with 95% CI.
# This is useful for statistical checking even if the manuscript retains the
# Cohen's d presentation.

fig2a_plot_data <- fig2a_overall %>%
  filter(Metric == "Beta_raw") %>%
  mutate(
    Group = factor(Group, levels = GROUP_ORDER),
    Fig2_treatment = factor(Fig2_treatment, levels = c("W", "H", "D", "C")),
    Star_label = ifelse(
      Significance_Holm_four_overall %in% c("*", "**", "***"),
      Significance_Holm_four_overall,
      ""
    )
  )

if (nrow(fig2a_plot_data) > 0) {
  # Dynamic star offset so negative and positive effects are both readable.
  star_offset <- max(abs(fig2a_plot_data$Cohen_d_overall), na.rm = TRUE) * 0.08
  if (!is.finite(star_offset) || star_offset == 0) star_offset <- 0.08

  fig2a_plot_data <- fig2a_plot_data %>%
    mutate(
      Star_y = Cohen_d_overall + ifelse(Cohen_d_overall >= 0, star_offset, -star_offset)
    )

  p_fig2a_cohen <- ggplot(
    fig2a_plot_data,
    aes(x = Fig2_treatment, y = Cohen_d_overall)
  ) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_col(width = 0.68) +
    geom_text(
      aes(y = Star_y, label = Star_label),
      size = 4.0
    ) +
    facet_wrap(~ Group, scales = "free_y", nrow = 1) +
    labs(
      x = NULL,
      y = "Overall Cohen's d"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold")
    )

  ggsave(
    file.path(FIGDIR, "Fig2a_OVERALL_Cohen_d_with_LMM_Holm_Beta_raw.pdf"),
    p_fig2a_cohen, width = 11, height = 3.8
  )
  ggsave(
    file.path(FIGDIR, "Fig2a_OVERALL_Cohen_d_with_LMM_Holm_Beta_raw.png"),
    p_fig2a_cohen, width = 11, height = 3.8, dpi = 300
  )

  p_fig2a_lmm <- ggplot(
    fig2a_plot_data,
    aes(x = Fig2_treatment, y = Estimate)
  ) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_errorbar(
      aes(ymin = CI_low, ymax = CI_high),
      width = 0.15
    ) +
    geom_point(aes(shape = Significance_Holm_four_overall), size = 2.6) +
    facet_wrap(~ Group, scales = "free_y", nrow = 1) +
    labs(
      x = NULL,
      y = expression("Overall LMM contrast " * Delta * beta),
      shape = "Holm significance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave(
    file.path(FIGDIR, "Fig2a_OVERALL_LMM_contrasts_Beta_raw.pdf"),
    p_fig2a_lmm, width = 11, height = 3.8
  )
  ggsave(
    file.path(FIGDIR, "Fig2a_OVERALL_LMM_contrasts_Beta_raw.png"),
    p_fig2a_lmm, width = 11, height = 3.8, dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 16. PLOTTING HELPERS
# ------------------------------------------------------------------------------

theme_revision <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.title = element_blank()
    )
}

save_pdf_png <- function(plot, filename, width, height) {
  
  ggsave(
    file.path(FIGDIR, paste0(filename, ".pdf")),
    plot = plot,
    width = width,
    height = height
  )
  
  ggsave(
    file.path(FIGDIR, paste0(filename, ".png")),
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 17. FIGURE A: 12 TREATMENTS x 4 BLOCKS FOR BETA_RAW
# ------------------------------------------------------------------------------

if ("Beta_raw" %in% names(dat)) {
  
  summ_beta <- dat %>%
    filter(is.finite(Beta_raw)) %>%
    group_by(Group, Treatment) %>%
    summarise(
      N = n(),
      Mean = mean(Beta_raw),
      SD = sd(Beta_raw),
      SE = SD / sqrt(N),
      .groups = "drop"
    )
  
  p_treat <- ggplot(
    dat %>% filter(is.finite(Beta_raw)),
    aes(
      x = Treatment,
      y = Beta_raw
    )
  ) +
    geom_point(
      aes(shape = block),
      position = position_jitter(
        width = 0.08,
        height = 0
      ),
      size = 2,
      alpha = 0.75
    ) +
    geom_point(
      data = summ_beta,
      aes(
        x = Treatment,
        y = Mean
      ),
      inherit.aes = FALSE,
      size = 2.7
    ) +
    geom_errorbar(
      data = summ_beta,
      aes(
        x = Treatment,
        ymin = Mean - SE,
        ymax = Mean + SE
      ),
      inherit.aes = FALSE,
      width = 0.15
    ) +
    facet_wrap(
      ~ Group,
      scales = "free_y",
      nrow = 1
    ) +
    labs(
      x = NULL,
      y = expression(beta[t])
    ) +
    theme_revision() +
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )
  
  save_pdf_png(
    p_treat,
    "Fig_QC_12treatments_Beta_raw",
    width = 12,
    height = 3.4
  )
}

# ------------------------------------------------------------------------------
# 18. FIGURE S13a: RAW VS DETRENDED BETA
# ------------------------------------------------------------------------------

if (
  all(
    c(
      "Beta_raw",
      "Beta_detrended"
    ) %in% names(dat)
  )
) {
  
  p_det <- ggplot(
    dat,
    aes(
      x = Beta_raw,
      y = Beta_detrended
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2
    ) +
    geom_point(
      aes(shape = block),
      size = 2,
      alpha = 0.75
    ) +
    facet_wrap(
      ~ Group,
      scales = "free",
      nrow = 1
    ) +
    labs(
      x = "Raw temporal TPL exponent",
      y = "Detrended temporal TPL exponent"
    ) +
    theme_revision()
  
  save_pdf_png(
    p_det,
    "Fig_S13a_raw_vs_detrended_beta",
    width = 11,
    height = 3.2
  )
}

# ------------------------------------------------------------------------------
# 19. FIGURE S13b: OVERALL WARMING EFFECT ON RAW / DETRENDED BETA
# ------------------------------------------------------------------------------

warming_beta <- primary_overall %>%
  filter(
    Family == "Warming",
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  mutate(
    Metric_label = dplyr::recode(
      as.character(Metric),
      Beta_raw = "Raw beta",
      Beta_detrended = "Detrended beta"
    ),
    Group = factor(
      Group,
      levels = GROUP_ORDER
    )
  )

if (nrow(warming_beta) > 0) {
  
  p_wbeta <- ggplot(
    warming_beta,
    aes(
      x = Estimate,
      y = Group,
      shape = Metric_label
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      width = 0.15,
      position = position_dodge(
        width = 0.5
      )
    ) +
    geom_point(
      size = 2.6,
      position = position_dodge(
        width = 0.5
      )
    ) +
    labs(
      x = "Overall warmed - unwarmed effect",
      y = NULL
    ) +
    theme_revision() +
    theme(
      legend.position = "top"
    )
  
  save_pdf_png(
    p_wbeta,
    "Fig_S13b_overall_warming_raw_detrended_beta",
    width = 6.5,
    height = 4.5
  )
}

# ------------------------------------------------------------------------------
# 20. FIGURE S13c: REFERENCE-ABUNDANCE VARIANCE
# ------------------------------------------------------------------------------

refvar_warm <- primary_overall %>%
  filter(
    Family == "Warming",
    Metric == "LogVar_at_ref_raw"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = GROUP_ORDER
    )
  )

if (nrow(refvar_warm) > 0) {
  
  p_ref <- ggplot(
    refvar_warm,
    aes(
      x = Estimate,
      y = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      width = 0.15
    ) +
    geom_point(size = 2.6) +
    labs(
      x = "Overall warming effect on log variance at reference abundance",
      y = NULL
    ) +
    theme_revision()
  
  save_pdf_png(
    p_ref,
    "Fig_S13c_reference_variance_overall_warming",
    width = 6.3,
    height = 4.2
  )
}

# ------------------------------------------------------------------------------
# 21. FIGURE S14: HELLINGER COMMUNITY TEMPORAL VARIABILITY
# ------------------------------------------------------------------------------

ctv_overall <- primary_overall %>%
  filter(
    Family == "Warming",
    Metric == "CTV_Hellinger_RMS"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = GROUP_ORDER
    )
  )

if (nrow(ctv_overall) > 0) {
  
  p_ctv <- ggplot(
    ctv_overall,
    aes(
      x = Estimate,
      y = Group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      width = 0.15
    ) +
    geom_point(size = 2.6) +
    labs(
      x = "Overall warming effect on Hellinger temporal dispersion",
      y = NULL
    ) +
    theme_revision()
  
  save_pdf_png(
    p_ctv,
    "Fig_S14a_CTV_overall_warming",
    width = 6.3,
    height = 4.2
  )
}

ctv_context <- primary_contrasts %>%
  filter(
    Family == "Warming",
    Metric == "CTV_Hellinger_RMS"
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = GROUP_ORDER
    ),
    Context = factor(
      Context,
      levels = rev(
        contrast_def %>%
          filter(Family == "Warming") %>%
          arrange(Contrast_order) %>%
          pull(Context)
      )
    )
  )

if (nrow(ctv_context) > 0) {
  
  p_ctv2 <- ggplot(
    ctv_context,
    aes(
      x = Estimate,
      y = Context
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      width = 0.15
    ) +
    geom_point(
      aes(
        shape = Significance_Holm
      ),
      size = 2.2
    ) +
    facet_wrap(
      ~ Group,
      scales = "free_x",
      nrow = 1
    ) +
    labs(
      x = "Warmed - unwarmed effect on Hellinger temporal dispersion",
      y = NULL
    ) +
    theme_revision() +
    theme(
      legend.position = "top"
    )
  
  save_pdf_png(
    p_ctv2,
    "Fig_S14b_CTV_context_specific_warming",
    width = 13,
    height = 4.8
  )
}

# ------------------------------------------------------------------------------
# 22. EFFECT-SPECIFIC BETA_RAW FOREST PLOTS
# ------------------------------------------------------------------------------

if (
  "Beta_raw" %in% unique(
    primary_contrasts$Metric
  )
) {
  
  for (fam in unique(contrast_def$Family)) {
    
    pd <- primary_contrasts %>%
      filter(
        Metric == "Beta_raw",
        Family == fam
      ) %>%
      mutate(
        Group = factor(
          Group,
          levels = GROUP_ORDER
        ),
        Context = factor(
          Context,
          levels = rev(
            contrast_def %>%
              filter(Family == fam) %>%
              arrange(Contrast_order) %>%
              pull(Context)
          )
        )
      )
    
    if (nrow(pd) == 0) next
    
    p <- ggplot(
      pd,
      aes(
        x = Estimate,
        y = Context
      )
    ) +
      geom_vline(
        xintercept = 0,
        linetype = 2
      ) +
      geom_errorbar(
        aes(
          xmin = CI_low,
          xmax = CI_high
        ),
        width = 0.14
      ) +
      geom_point(
        aes(
          shape = Significance_Holm
        ),
        size = 2.2
      ) +
      facet_wrap(
        ~ Group,
        scales = "free_x",
        nrow = 1
      ) +
      labs(
        x = paste0(
          fam,
          " effect on temporal TPL beta"
        ),
        y = NULL
      ) +
      theme_revision() +
      theme(
        legend.position = "top"
      )
    
    safe_fam <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      fam
    )
    
    save_pdf_png(
      p,
      paste0(
        "Fig_effect_specific_Beta_raw_",
        safe_fam
      ),
      width = 13,
      height = ifelse(
        fam %in% c(
          "Half_precipitation",
          "Double_precipitation"
        ),
        4.2,
        4.8
      )
    )
  }
}

# ------------------------------------------------------------------------------
# 23. COMPACT TABLE FOR REVIEWER RESPONSE
# ------------------------------------------------------------------------------

reviewer_compact <- unified_overall %>%
  filter(
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended",
      "LogVar_at_ref_raw",
      "CTV_Hellinger_RMS"
    )
  ) %>%
  arrange(
    Metric,
    Family,
    Group
  )

write_out(
  reviewer_compact,
  "14_REVIEWER_compact_overall_effects.csv"
)

# For Beta_raw, keep all planned context-specific effects with all evidence.
reviewer_beta_context <- unified %>%
  filter(
    Metric %in% c(
      "Beta_raw",
      "Beta_detrended"
    )
  ) %>%
  arrange(
    Metric,
    Family,
    Group,
    Contrast_order
  )

write_out(
  reviewer_beta_context,
  "14_REVIEWER_beta_effect_specific_all_evidence.csv"
)

# ------------------------------------------------------------------------------
# 24. README / INTERPRETATION GUIDE
# ------------------------------------------------------------------------------

readme <- c(
  "Reviewer 2 revised analysis: categorical precipitation primary model",
  "",
  "DATA STRUCTURE",
  "For each microbial group: 48 subplot-level metric values = 12 treatments x 4 blocks.",
  "There are 24 physical whole plots nested within 4 blocks, with clipped/unclipped subplots within each whole plot.",
  "",
  "PRIMARY MODEL",
  "metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)",
  "WarmF and ClipF are two-level factors; PrecipF is a three-level factor (Ambient/Half/Double).",
  "This primary model estimates Half-vs-Ambient and Double-vs-Ambient effects independently.",
  "",
  "CONTINUOUS-GRADIENT SENSITIVITY",
  "metric ~ Warming * Precipitation * Clipping + (1 | Block/WholePlotID)",
  "Precipitation is coded -0.5/0/1, corresponding to multipliers 0.5/1/2 centered at ambient.",
  "This model imposes a linear precipitation-response restriction and is NOT the primary model.",
  "03_MODEL_COMPARISON_categorical_vs_continuous.csv reports ML AIC/BIC/logLik/LRT comparisons.",
  "03_MODEL_COMPARISON_leave_one_block_out*.csv reports optional block-wise predictive comparison for beta metrics.",
  "",
  "MODEL QUALITY / DIAGNOSTICS",
  "03_PRIMARY_categorical_model_diagnostics.csv reports singular fits, variance components, marginal/conditional R2, residual summaries, QQ correlation, and an absolute-residual-vs-fitted diagnostic P value.",
  "Residual-vs-fitted and QQ plots are saved for Beta_raw/Beta_detrended under figures/Diag_*.",
  "R2 is descriptive and is not a pass/fail threshold for validity of randomized treatment contrasts.",
  "The critical checks are correct experimental hierarchy, reasonable residual behavior, non-pathological singularity, and robustness to the block-fixed sensitivity.",
  "",
  "PRIMARY CONTRASTS",
  "04_PRIMARY_effect_specific_contrasts_Holm.csv contains all 20 biological context-specific contrasts from the categorical primary model.",
  "Warming: W-Ctrl, WH-H, WD-D, WC-C, WHC-HC, WDC-DC",
  "Clipping: C-Ctrl, WC-W, HC-H, DC-D, WHC-WH, WDC-WD",
  "Half precipitation: H-Ctrl, WH-W, HC-C, WHC-WC",
  "Double precipitation: D-Ctrl, WD-W, DC-C, WDC-WC",
  "Holm corrections are provided within each biological family and across all 20 contrasts.",
  "",
  "FIGURE 2a OVERALL W/H/D/C FACTOR EFFECTS",
  "Fig. 2a v3 uses overall factor effects, not W/H/D/C single treatment cells versus Ctrl.",
  "07_OVERALL_Cohen_d_Hedges_g.csv contains descriptive marginal pooled-SD Cohen d and Hedges g for overall Warming, Half, Double, and Clipping.",
  "05_PRIMARY_overall_directional_effects.csv contains the matching equal-weight categorical-LMM overall contrasts.",
  "P_Holm_four_overall adjusts the four overall effects within each Group x Metric.",
  "08_FIG2a_OVERALL_W_H_D_C_effectsize_and_LMM.csv merges overall Cohen d/Hedges g with LMM estimate, CI, raw P, and Holm P.",
  "08_FIG2a_OVERALL_beta_compact.csv is the compact Beta_raw/Beta_detrended version for manuscript checking.",
  "figures/Fig2a_OVERALL_Cohen_d_with_LMM_Holm_Beta_raw.* is the primary Cohen-d style Fig. 2a candidate.",
  "figures/Fig2a_OVERALL_LMM_contrasts_Beta_raw.* is an optional LMM Delta-beta + 95% CI check plot.",
  "",
  "SECONDARY SINGLE-TREATMENT DIAGNOSTIC",
  "04_PRIMARY_single_treatment_vs_Ctrl_Holm.csv retains W-Ctrl, H-Ctrl, D-Ctrl, C-Ctrl for transparency only.",
  "07_single_treatment_Cohen_d_vs_LMM_contrast_diagnostic.csv compares those raw cell differences with categorical-LMM cell contrasts.",
  "",
  "OVERALL EFFECT DEFINITIONS",
  "Warming = equal-weight mean of W-Ctrl, WH-H, WD-D, WC-C, WHC-HC, WDC-DC.",
  "Half = equal-weight mean of H-Ctrl, WH-W, HC-C, WHC-WC.",
  "Double = equal-weight mean of D-Ctrl, WD-W, DC-C, WDC-WC.",
  "Clipping = equal-weight mean of C-Ctrl, WC-W, HC-H, DC-D, WHC-WH, WDC-WD.",
  "",
  "CONTINUOUS SENSITIVITY OUTPUTS",
  "05_CONTINUOUS_sensitivity_effect_specific_contrasts.csv reports the same treatment contrasts under the linear precipitation-gradient restriction.",
  "05_MODEL_COMPARISON_effect_specific_categorical_vs_continuous.csv compares categorical and continuous contrast estimates directly.",
  "05_CONTINUOUS_sensitivity_precipitation_slopes.csv reports precipitation slopes from the numeric model.",
  "",
  "EFFECT SIZE",
  "07_Cohen_d_Hedges_g_effect_specific.csv contains context-specific descriptive pooled-SD Cohen d and Hedges g.",
  "07_OVERALL_Cohen_d_Hedges_g.csv contains the overall marginal descriptive Cohen d and Hedges g used for Fig. 2a.",
  "These standardized effect sizes are not used as the basis for significance; stars/P values come from the categorical hierarchical LMM.",
  "",
  "BETA-SE PROPAGATION",
  "10_Beta_SE_propagation_* files propagate SE_Beta^2 using the same categorical fixed-effect structure as the primary model.",
  "This is a sensitivity analysis and does not replace the primary LMM.",
  "",
  "BETA-COMMUNITY VARIABILITY ASSOCIATION",
  "13_beta_CTV_association_adjusted_descriptive.csv uses:",
  "CTV_Hellinger_RMS ~ Beta_raw_z + WarmF * PrecipF * ClipF + (1 | Block/WholePlotID).",
  "",
  "BETA / REFERENCE-VARIANCE CO-VARIATION UNDER WARMING",
  "The adjusted change-change model uses categorical precipitation/clipping covariates:",
  "Delta_LogVar_ref ~ Delta_Beta + PrecipF * ClipF + Block + (1 | PairID).",
  "",
  "REPORTING PRIORITY",
  "1) Categorical hierarchical mixed model and its overall/context-specific planned contrasts.",
  "2) Fig. 2a: overall Cohen d/Hedges g for descriptive magnitude + matching overall LMM Holm P for inference.",
  "3) Model-quality diagnostics and categorical-vs-continuous sensitivity comparison.",
  "4) Context-specific Cohen d / Hedges g as descriptive magnitude/direction only.",
  "5) Block-matched and response-ratio analyses as sensitivity.",
  "6) Beta-SE multilevel meta-regression as uncertainty-propagation sensitivity."
)

writeLines(
  readme,
  file.path(OUTDIR, "README_interpretation.txt")
)

# ------------------------------------------------------------------------------
# 25. FINAL CONSOLE SUMMARY
# ------------------------------------------------------------------------------

cat("\n")
cat("====================================================================\n")
cat("Reviewer 2 revised categorical-primary analysis completed\n")
cat("====================================================================\n\n")
cat("Input:\n  ", INPUT_FILE, "\n\n", sep = "")
cat("Output:\n  ", normalizePath(OUTDIR, mustWork = FALSE), "\n\n", sep = "")

cat("Model comparison for Beta_raw/Beta_detrended:\n")
print(
  model_comparison_out %>%
    filter(Metric %in% c("Beta_raw", "Beta_detrended")) %>%
    select(
      Group, Metric,
      Delta_AIC_cat_minus_cont,
      Delta_BIC_cat_minus_cont,
      LR_Chisq, LR_df, LR_P,
      Preferred_by_AIC
    )
)

cat("\nPrimary categorical model diagnostics for beta metrics:\n")
print(
  diag_out %>%
    filter(Metric %in% c("Beta_raw", "Beta_detrended")) %>%
    select(
      Group, Metric, Singular,
      R2_marginal, R2_conditional,
      Max_abs_standardized_residual,
      QQ_correlation,
      AbsResidual_vs_Fitted_P
    )
)

cat("\nSecondary W/H/D/C-vs-Ctrl diagnostic contrasts for Beta_raw:\n")
print(
  fig2_single4 %>%
    filter(Metric == "Beta_raw") %>%
    select(
      Group, Fig2_treatment,
      Estimate, SE, CI_low, CI_high,
      P_raw, P_Holm_single4,
      Significance_Holm_single4
    )
)

cat("\nPRIMARY Fig. 2a overall effects for Beta_raw:\n")
print(
  fig2a_overall %>%
    filter(Metric == "Beta_raw") %>%
    select(
      Group, Fig2_treatment,
      Cohen_d_overall, Hedges_g_overall,
      Estimate, SE, CI_low, CI_high,
      P_raw, P_Holm_four_overall,
      Significance_Holm_four_overall
    )
)

if (nrow(meta_overall) > 0) {
  cat("\nBeta-SE propagated overall effects:\n")
  print(
    meta_overall %>%
      select(
        Group, Metric, Family,
        Estimate_SEprop, SE_SEprop,
        P_SEprop, P_SEprop_Holm_four_overall
      )
  )
}

cat("\nIMPORTANT:\n")
cat("  The categorical precipitation model is primary because H and D are\n")
cat("  interpreted as distinct experimental treatments. The continuous\n")
cat("  precipitation model is retained as a constrained sensitivity analysis.\n")
cat("  Fig. 2a v3 uses OVERALL W/H/D/C factor effects: descriptive overall\n")
cat("  Cohen's d/Hedges' g for magnitude and matching categorical-LMM overall\n")
cat("  contrasts (Holm-adjusted across four effects) for formal inference.\n")
cat("  Inspect model-comparison and diagnostic outputs before final manuscript wording.\n\n")

# Restore user's contrast options.
options(old_contrasts)

# ==============================================================================
# tTTPL ADD-ON: ALL 11 TREATMENTS VS CTRL FOR TABLE 1
#
# Paste/run this block AFTER the primary categorical models have been fitted in:
#   3.Reviewer2_complete_supplementary_statistics_v3_overall_Fig2.R
#
# It expects the following objects/functions from that script to already exist:
#   dat, primary_models, TREATMENT_LEVELS, GROUP_ORDER, MAIN_METRICS,
#   PRIMARY_MODEL_LABEL, make_emm_weight, label_factor_emm_treatment,
#   cohen_pooled, sig_label, write_out
#
# Primary inference:
#   metric ~ WarmF * PrecipF * ClipF + (1 | Block/WholePlotID)
#
# Outputs:
#   06_PRIMARY_all_treatments_vs_Ctrl.csv
#   07_Cohen_d_Hedges_g_vs_Ctrl.csv
#   08_TABLE1_tTTPL_Beta_raw_vs_Ctrl.csv
#
# Multiplicity:
#   Holm correction across the 11 treatment-vs-Ctrl contrasts separately
#   within each Group x Metric.
# ============================================================================== 

# ------------------------------------------------------------------------------
# A. CATEGORICAL-LMM CONTRASTS: ALL 11 TREATMENTS VS CTRL
# ------------------------------------------------------------------------------

extract_all_treatments_vs_ctrl <- function(model, group, metric) {
  
  emm <- emmeans::emmeans(
    model,
    specs = ~ WarmF * PrecipF * ClipF,
    lmer.df = "satterthwaite"
  )
  
  emm_grid <- as.data.frame(emm)
  
  emm_levels <- label_factor_emm_treatment(
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
      "Could not reconstruct the unique 12-treatment categorical EMM grid for ",
      group, " / ", metric, "."
    )
  }
  
  rows <- list()
  k <- 1
  
  for (trt in setdiff(TREATMENT_LEVELS, "Ctrl")) {
    
    w <- make_emm_weight(
      emm_levels = emm_levels,
      target = trt,
      reference = "Ctrl"
    )
    
    contrast_label <- paste0(trt, " - Ctrl")
    
    cc <- emmeans::contrast(
      emm,
      method = setNames(
        list(w),
        contrast_label
      ),
      adjust = "none"
    )
    
    ss <- as.data.frame(
      summary(
        cc,
        infer = c(TRUE, TRUE),
        adjust = "none"
      )
    )
    
    rows[[k]] <- tibble::tibble(
      Group = group,
      Metric = metric,
      Treatment = trt,
      Treatment_order = match(trt, TREATMENT_LEVELS),
      Reference = "Ctrl",
      Contrast = contrast_label,
      Model_basis = PRIMARY_MODEL_LABEL,
      Estimate = ss$estimate[1],
      SE = ss$SE[1],
      df = ss$df[1],
      CI_low = ss$lower.CL[1],
      CI_high = ss$upper.CL[1],
      t_ratio = ss$t.ratio[1],
      P_raw = ss$p.value[1]
    )
    
    k <- k + 1
  }
  
  bind_rows(rows) %>%
    group_by(Group, Metric) %>%
    mutate(
      # The Table-1 multiplicity family is the 11 treatment-vs-Ctrl contrasts.
      P_Holm_11_vsCtrl = p.adjust(
        P_raw,
        method = "holm"
      ),
      Significance_Holm_11_vsCtrl = sig_label(
        P_Holm_11_vsCtrl
      ),
      Significant_Holm_11_vsCtrl =
        !is.na(P_Holm_11_vsCtrl) & P_Holm_11_vsCtrl < 0.05
    ) %>%
    ungroup() %>%
    arrange(Group, Metric, Treatment_order)
}


all_vsctrl_rows <- list()
idx_vsctrl <- 1

for (key in names(primary_models)) {
  
  ans <- primary_models[[key]]
  
  if (is.null(ans$categorical)) {
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
    extract_all_treatments_vs_ctrl(
      model = ans$categorical,
      group = g,
      metric = metric
    ),
    error = function(e) e
  )
  
  if (!inherits(rr, "error")) {
    all_vsctrl_rows[[idx_vsctrl]] <- rr
    idx_vsctrl <- idx_vsctrl + 1
  } else {
    warning(
      "All-treatment-vs-Ctrl contrast extraction failed for ",
      key,
      ": ",
      conditionMessage(rr)
    )
  }
}

all_treatments_vs_ctrl <- bind_rows(
  all_vsctrl_rows
)

write_out(
  all_treatments_vs_ctrl,
  "06_PRIMARY_all_treatments_vs_Ctrl.csv"
)


# ------------------------------------------------------------------------------
# B. DESCRIPTIVE COHEN'S d / HEDGES' g: ALL 11 TREATMENTS VS CTRL
# ------------------------------------------------------------------------------

cohen_vsctrl_rows <- list()
idx_cohen_vsctrl <- 1

for (g in unique(as.character(na.omit(dat$Group)))) {
  
  dg <- dat %>%
    filter(Group == g)
  
  for (metric in MAIN_METRICS) {
    
    if (!metric %in% names(dg)) {
      next
    }
    
    x_ctrl <- dg %>%
      filter(Treatment == "Ctrl") %>%
      pull(all_of(metric))
    
    for (trt in setdiff(TREATMENT_LEVELS, "Ctrl")) {
      
      x_trt <- dg %>%
        filter(Treatment == trt) %>%
        pull(all_of(metric))
      
      es <- cohen_pooled(
        x = x_trt,
        y = x_ctrl
      )
      
      cohen_vsctrl_rows[[idx_cohen_vsctrl]] <- tibble::tibble(
        Group = g,
        Metric = metric,
        Treatment = trt,
        Treatment_order = match(trt, TREATMENT_LEVELS),
        Reference = "Ctrl",
        Contrast = paste0(trt, " - Ctrl"),
        N_target = unname(es["n_target"]),
        N_reference = unname(es["n_reference"]),
        Mean_target = unname(es["mean_target"]),
        Mean_reference = unname(es["mean_reference"]),
        Raw_difference =
          unname(es["mean_target"]) - unname(es["mean_reference"]),
        SD_target = unname(es["sd_target"]),
        SD_reference = unname(es["sd_reference"]),
        SD_pooled = unname(es["sd_pooled"]),
        Cohen_d_pooled = unname(es["d"]),
        Hedges_g = unname(es["g"])
      )
      
      idx_cohen_vsctrl <- idx_cohen_vsctrl + 1
    }
  }
}

cohen_vsctrl <- bind_rows(
  cohen_vsctrl_rows
) %>%
  arrange(Group, Metric, Treatment_order)

write_out(
  cohen_vsctrl,
  "07_Cohen_d_Hedges_g_vs_Ctrl.csv"
)


# ------------------------------------------------------------------------------
# C. TABLE 1 BASE: BETA_raw ONLY
#
# Numeric value shown in Table 1:
#   Cohen_d_pooled
#
# Bold formatting later:
#   Significant_Holm_11_vsCtrl == TRUE
#
# LMM Estimate/CI/P are retained here for auditability but do not need to be
# printed in the final compact manuscript table.
# ------------------------------------------------------------------------------

table1_tttpl_beta_raw <- cohen_vsctrl %>%
  filter(Metric == "Beta_raw") %>%
  select(
    Group,
    Metric,
    Treatment,
    Treatment_order,
    Reference,
    Contrast,
    N_target,
    N_reference,
    Mean_target,
    Mean_reference,
    Raw_difference,
    Cohen_d_pooled,
    Hedges_g
  ) %>%
  left_join(
    all_treatments_vs_ctrl %>%
      filter(Metric == "Beta_raw") %>%
      select(
        Group,
        Metric,
        Treatment,
        Estimate,
        SE,
        df,
        CI_low,
        CI_high,
        t_ratio,
        P_raw,
        P_Holm_11_vsCtrl,
        Significance_Holm_11_vsCtrl,
        Significant_Holm_11_vsCtrl
      ),
    by = c(
      "Group",
      "Metric",
      "Treatment"
    )
  ) %>%
  mutate(
    # In a balanced complete categorical factorial model the LMM point estimate
    # should be identical, or numerically almost identical, to the raw cell-mean
    # difference. This column is a useful QC check.
    LMM_minus_raw_difference = Estimate - Raw_difference,
    Same_direction_raw_and_LMM = case_when(
      is.na(Estimate) | is.na(Raw_difference) ~ NA,
      Estimate == 0 | Raw_difference == 0 ~ TRUE,
      TRUE ~ sign(Estimate) == sign(Raw_difference)
    )
  ) %>%
  arrange(
    factor(Group, levels = GROUP_ORDER),
    Treatment_order
  )

write_out(
  table1_tttpl_beta_raw,
  "08_TABLE1_tTTPL_Beta_raw_vs_Ctrl.csv"
)


# ------------------------------------------------------------------------------
# D. QC CHECKS
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("TABLE 1 tTTPL: 11 TREATMENTS VS CTRL\n")
cat("============================================================\n\n")

cat(
  "Expected Beta_raw rows = ",
  length(unique(as.character(na.omit(dat$Group)))) *
    length(setdiff(TREATMENT_LEVELS, "Ctrl")),
  "\n",
  sep = ""
)

cat(
  "Observed Beta_raw rows = ",
  nrow(
    table1_tttpl_beta_raw
  ),
  "\n\n",
  sep = ""
)

print(
  table1_tttpl_beta_raw %>%
    select(
      Group,
      Treatment,
      Cohen_d_pooled,
      Estimate,
      P_raw,
      P_Holm_11_vsCtrl,
      Significance_Holm_11_vsCtrl,
      Significant_Holm_11_vsCtrl,
      LMM_minus_raw_difference
    ),
  n = Inf
)

if (
  any(
    abs(table1_tttpl_beta_raw$LMM_minus_raw_difference) > 1e-8,
    na.rm = TRUE
  )
) {
  warning(
    "At least one categorical-LMM treatment-vs-Ctrl point estimate differs ",
    "from the raw cell-mean difference by more than 1e-8. Inspect the Table-1 QC output."
  )
}

cat("\nFiles written:\n")
cat("  06_PRIMARY_all_treatments_vs_Ctrl.csv\n")
cat("  07_Cohen_d_Hedges_g_vs_Ctrl.csv\n")
cat("  08_TABLE1_tTTPL_Beta_raw_vs_Ctrl.csv\n")
