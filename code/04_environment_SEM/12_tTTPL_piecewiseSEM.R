# ==============================================================================
# Fig. 3b-c: taxonomic temporal TPL (tTTPL)
# Reduced mechanism-led mixed-effects piecewise SEM respecting the split-plot
# design. The environmental backbone and beta equations are fixed before
# fitting; no forward/backward/stepwise selection occurs in the engine.
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

DATA_FILE <- "./results/10_environment_SEM_inputs/tTTPL_SEM_input.csv"
META_FILE <- "./results/10_environment_SEM_inputs/treatment_metadata_48plots.csv"
OUT_DIR <- "./results/11_tTTPL_piecewiseSEM"

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
# Archaea additionally contains the correlated error:
#   TC_mean %~~% TC_cv
# This term is handled by the engine as a piecewiseSEM formula.cerror and is
# used in the d-separation/Fisher's C calculation, but not multiplied into
# directed indirect/total effects.
MODEL_SPEC <- make_tTTPL_model_spec()

N_UNIQUE_FIXED_EDGES <- nrow(unique(do.call(rbind, lapply(
  MODEL_SPEC, function(x) spec_edges(x$equations)
))[, c("From", "To"), drop = FALSE]))

TARGET_LABELS <- c(
  Bac_beta = "Bacteria",
  Arch_beta = "Archaea",
  Fun_18S_beta = "Fungi",
  Pro_beta = "Protists",
  Nem_beta = "Nematoda"
)

# ------------------------------- RUN ------------------------------------------
result_tTTPL <- run_splitplot_piecewise_sem(
  data_file = DATA_FILE,
  meta_file = META_FILE,
  out_dir = OUT_DIR,
  analysis_label = "tTTPL",
  figure_label = "Fig3",
  model_spec = MODEL_SPEC,
  target_labels = TARGET_LABELS,
  model_spec_basis = paste(
    "reduced mechanism-led environmental backbone shared with pTTPL;",
    "taxon-specific tTTPL beta equations informed by prior LMM and reduced-SEM diagnostics;",
    "TC_mean ~~ TC_cv specified as an Archaea correlated error;",
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
