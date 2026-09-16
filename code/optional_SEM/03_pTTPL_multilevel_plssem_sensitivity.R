# ==============================================================================
# Supplementary Fig. 11b-c pTTPL: optional multilevel plssem sensitivity analysis
# Primary hierarchy is (1|block/WholePlotID). Current plssem accepts one random
# cluster only, so this sensitivity model fits block fixed contrasts plus
# (1|BlockWholePlotID). See README_CN.md; this is not the primary model.
# Do not use as the primary manuscript analysis without additional validation.
# ==============================================================================

rm(list = ls())
# ------------------------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------------------------
# Run this script from the repository root, or set TPL_PROJECT_ROOT to the
# project directory that contains the data folders documented in README.md.
PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
source(file.path("code", "04_environment_SEM", "11_splitplot_piecewiseSEM_engine.R"))
source(file.path("code", "optional_SEM", "01_splitplot_plssem_engine.R"))
source(file.path("code", "04_environment_SEM", "10_TPL_SEM_model_specs.R"))

DATA_FILE <- "./results/10_environment_SEM_inputs/pTTPL_SEM_input.csv"
META_FILE <- "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"
OUT_DIR <- "./results/optional_pTTPL_plssem_sensitivity"

RUN_BOOTSTRAP <- TRUE
N_BOOT <- 999L
RUN_LOBO <- TRUE
SEED <- 8302026L
N_CORES <- if (.Platform$OS.type == "windows") 1L else {
  max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
}

MODEL_SPEC <- make_pTTPL_model_spec()
TARGET_LABELS <- c(
  Bac_pTTPL_beta = "Bacteria", Arch_pTTPL_beta = "Archaea",
  Fun_18S_pTTPL_beta = "Fungi", Pro_pTTPL_beta = "Protists",
  Nem_pTTPL_beta = "Nematoda"
)

result_pTTPL_plssem <- run_splitplot_plssem(
  data_file = DATA_FILE,
  meta_file = META_FILE,
  out_dir = OUT_DIR,
  analysis_label = "pTTPL",
  figure_label = "SuppFig11_plssem_sensitivity",
  model_spec = MODEL_SPEC,
  target_labels = TARGET_LABELS,
  run_bootstrap = RUN_BOOTSTRAP,
  n_boot = N_BOOT,
  seed = SEED,
  n_cores = N_CORES,
  run_lobo = RUN_LOBO
)
