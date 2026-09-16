# ==============================================================================
# Supplementary Fig. 11b-c: phylogenetic temporal TPL (pTTPL)
# Mixed-effects piecewise SEM respecting the split-plot design.
#
# The upstream environmental backbone is exactly the same as the current tTTPL
# model. Only the taxon-specific *_pTTPL_beta equations differ and are intended
# for diagnostic refinement after this first reduced pTTPL run.
# ==============================================================================

rm(list = ls())

# ------------------------------- SETTINGS -------------------------------------
# ------------------------------------------------------------------------------
# PROJECT ROOT
# ------------------------------------------------------------------------------
# Run this script from the repository root, or set TPL_PROJECT_ROOT to the
# project directory that contains the data folders documented in README.md.
PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)
SCRIPT_DIR <- file.path("code", "04_environment_SEM")
source(file.path(SCRIPT_DIR, "11_splitplot_piecewiseSEM_engine.R"))
source(file.path(SCRIPT_DIR, "10_TPL_SEM_model_specs.R"))

DATA_FILE <- "./results/10_environment_SEM_inputs/pTTPL_SEM_input.csv"
META_FILE <- "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"
OUT_DIR <- "./results/12_pTTPL_piecewiseSEM"

RUN_BOOTSTRAP <- TRUE
N_BOOT <- 999L                 # use 99 for a quick test; 999/4999 for final
RUN_LOBO <- TRUE               # robustness only; never selects or deletes paths
SEED <- 8302026L
N_CORES <- if (.Platform$OS.type == "windows") 1L else {
  max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))
}
ALPHA_DISPLAY <- 0.05

TREATMENTS <- c("warm", "drying", "wetting", "clip")

# Consistent ordering for Fig. c / Fig. S9c source-effect comparisons.
# The engine retains every treatment and every intermediate/environmental
# variable in the clean *_all_source_effects_for_plot.csv output.
PLOT_SOURCE_ORDER <- c(
  "warm", "drying", "wetting", "clip",
  "temperature_mean",
  "moisture_mean", "moisture_cv",
  "pH_mean",
  "NO3N_cv",
  "TC_mean", "TC_cv",
  "plant_biomass_mean",
  "plant_richness_mean"
)

# ------------------------- FIXED MODEL SPECIFICATION ---------------------------
# Upstream equations come from make_shared_TPL_environmental_equations() and are
# therefore identical to tTTPL. The current pTTPL beta equations are reduced
# candidates and should be evaluated after this run using:
# Fisher's C, d-separation claims, target R2, whole-plot bootstrap and LOBO.
MODEL_SPEC <- make_pTTPL_model_spec()

N_UNIQUE_FIXED_EDGES <- nrow(unique(do.call(rbind, lapply(
  MODEL_SPEC, function(x) spec_edges(x$equations)
))[, c("From", "To"), drop = FALSE]))

TARGET_LABELS <- c(
  Bac_pTTPL_beta = "Bacteria",
  Arch_pTTPL_beta = "Archaea",
  Fun_18S_pTTPL_beta = "Fungi",
  Pro_pTTPL_beta = "Protists",
  Nem_pTTPL_beta = "Nematoda"
)

# ------------------------------- RUN ------------------------------------------
result_pTTPL <- run_splitplot_piecewise_sem(
  data_file = DATA_FILE,
  meta_file = META_FILE,
  out_dir = OUT_DIR,
  analysis_label = "pTTPL",
  figure_label = "SuppFig11",
  model_spec = MODEL_SPEC,
  target_labels = TARGET_LABELS,
  model_spec_basis = paste(
    "environmental backbone fixed to the current tTTPL specification;",
    "reduced candidate taxon-specific pTTPL beta equations;",
    "NO3N represented by NO3N_cv in the shared core model;",
    "no path selection within piecewiseSEM"
  ),
  treatments = TREATMENTS,
  alpha_display = ALPHA_DISPLAY,
  run_bootstrap = RUN_BOOTSTRAP,
  n_boot = N_BOOT,
  seed = SEED,
  n_cores = N_CORES,
  run_lobo = RUN_LOBO,
  plot_source_order = PLOT_SOURCE_ORDER
)
