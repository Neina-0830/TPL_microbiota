# ==============================================================================
# Fig. 3b-c tTTPL: optional multilevel plssem sensitivity analysis
# Primary hierarchy is (1|block/WholePlotID). Current plssem accepts one random
# cluster only, so this sensitivity model fits block fixed contrasts plus
# (1|BlockWholePlotID). See README_CN.md; this is not the primary model.
# MODEL_SPEC contains all 106 edges inherited from tTTPL PLS forward selection.
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

DATA_FILE <- "./results/10_environment_SEM_inputs/tTTPL_SEM_input.csv"
META_FILE <- "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"
OUT_DIR <- "./results/optional_tTTPL_plssem_sensitivity"

RUN_BOOTSTRAP <- TRUE
N_BOOT <- 999L
RUN_LOBO <- TRUE
SEED <- 8302026L
N_CORES <- if (.Platform$OS.type == "windows") 1L else {
  max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
}

MODEL_SPEC <- make_tTTPL_model_spec()
N_UNIQUE_FIXED_EDGES <- nrow(unique(do.call(rbind, lapply(
  MODEL_SPEC, function(x) spec_edges(x$equations)
))[, c("From", "To"), drop = FALSE]))
stopifnot(N_UNIQUE_FIXED_EDGES == 106L)
TARGET_LABELS <- c(
  Bac_beta = "Bacteria", Arch_beta = "Archaea", Fun_18S_beta = "Fungi",
  Pro_beta = "Protists", Nem_beta = "Nematoda"
)

result_tTTPL_plssem <- run_splitplot_plssem(
  data_file = DATA_FILE,
  meta_file = META_FILE,
  out_dir = OUT_DIR,
  analysis_label = "tTTPL",
  figure_label = "Fig3_plssem_sensitivity",
  model_spec = MODEL_SPEC,
  target_labels = TARGET_LABELS,
  model_spec_basis = paste(
    "all 106 tTTPL edges inherited from PLS forward selection;",
    "post-selection/exploratory specification"
  ),
  run_bootstrap = RUN_BOOTSTRAP,
  n_boot = N_BOOT,
  seed = SEED,
  n_cores = N_CORES,
  run_lobo = RUN_LOBO
)
