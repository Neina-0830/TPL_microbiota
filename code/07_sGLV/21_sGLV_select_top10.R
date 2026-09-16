# ==============================================================================
# 21_sGLV_select_top10.R
# Combine TPL, community-composition and richness calibration discrepancies,
# reproduce the supplied group/TPL-specific weighting, and retain the ten
# lowest-loss parameter sets separately for tTTPL and pTTPL.
#
# Outputs:
#   Calibration/<TPL>/<Group>_all_parameter_loss.csv
#   Calibration/<TPL>/<Group>_top10_parameters.csv
#
# The top-10 tables are the direct inputs to 22_sGLV_process_intensity_simulate.py.
# ==============================================================================

rm(list = ls())

PROJECT_ROOT <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(Sys.getenv("TPL_SGLV_WORK_DIR",
                                 unset = file.path(PROJECT_ROOT, "results", "17_sGLV")),
                      winslash = "/", mustWork = FALSE)
DATA_SGLV <- normalizePath(Sys.getenv("TPL_SGLV_DATA_DIR",
                                      unset = file.path(PROJECT_ROOT, "data", "sGLV")),
                           winslash = "/", mustWork = FALSE)
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

SCRIPT_DIR <- file.path(PROJECT_ROOT, "code", "07_sGLV")
cfg <- read.csv(file.path(SCRIPT_DIR, "sGLV_group_config.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
weights <- read.csv(file.path(SCRIPT_DIR, "sGLV_calibration_weights.csv"),
                    stringsAsFactors = FALSE, check.names = FALSE)

groups_env <- Sys.getenv("TPL_SGLV_GROUPS", unset = "all")
GROUPS <- if (tolower(groups_env) == "all") cfg$Group else trimws(strsplit(groups_env, ",")[[1]])

param_cols <- c("Immigration", "Interaction", "Growth", "Noise")
seed_cols <- c(param_cols, "Seed")

read_empirical_beta <- function(group, tpl) {
  # Use the canonical plot-level empirical beta outputs produced earlier in this
  # repository, rather than legacy warming-comparison files.
  manuscript_group <- if (group == "Fungi_18S") "Fungi" else group
  f <- if (tpl == "tTTPL") {
    file.path(PROJECT_ROOT, "results", "01_tTTPL",
              "UPARSE_all_groups_plot_level_revision_metrics.csv")
  } else {
    file.path(PROJECT_ROOT, "results", "04_pTTPL",
              "pTTPL_all_groups_plot_level_full_precision.csv")
  }
  if (!file.exists(f)) stop("Missing empirical beta table: ", f)
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("Group", "Beta_raw", "warm") %in% names(x))) {
    stop("Empirical beta table lacks Group/Beta_raw/warm: ", f)
  }
  z <- x[x$Group == manuscript_group & is.finite(x$Beta_raw), , drop = FALSE]
  if (!nrow(z)) stop("No empirical beta rows for ", manuscript_group, " in ", f)
  c(
    warming = mean(z$Beta_raw[as.numeric(z$warm) == 1], na.rm = TRUE),
    unwarming = mean(z$Beta_raw[as.numeric(z$warm) == 0], na.rm = TRUE)
  )
}

read_sim_beta <- function(group, condition, tpl) {
  prefix <- if (tpl == "tTTPL") "" else "MPD_weighted_"
  f <- file.path(
    ROOT, "Taylor", group, "raw", condition,
    paste0(prefix,
      "SLM_all_immigration_interaction_growth_all_seed_noise_all_species_Taylor_beta.csv")
  )
  if (!file.exists(f)) stop("Missing calibration TPL table: ", f)
  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  x$Seed <- as.character(x$Seed)
  x
}

for (tpl in c("tTTPL", "pTTPL")) {
  for (group in GROUPS) {
    emp <- read_empirical_beta(group, tpl)
    w <- read_sim_beta(group, "warming", tpl)
    u <- read_sim_beta(group, "unwarming", tpl)

    pair <- merge(
      w[, c(seed_cols, "Beta")],
      u[, c(seed_cols, "Beta")],
      by = seed_cols,
      suffixes = c("_warming", "_unwarming")
    )
    if (!nrow(pair)) stop(group, " ", tpl, ": no matched warmed/unwarmed beta rows.")

    pair$Beta_diff <- (pair$Beta_warming - emp["warming"])^2 +
      (pair$Beta_unwarming - emp["unwarming"])^2
    pair$Delta_beta_diff <- (
      (pair$Beta_warming - pair$Beta_unwarming) -
      (emp["warming"] - emp["unwarming"])
    )^2

    beta_param <- aggregate(
      cbind(Beta_diff, Delta_beta_diff) ~ Immigration + Interaction + Growth + Noise,
      data = pair, FUN = mean, na.rm = TRUE
    )
    names(beta_param)[names(beta_param) == "Beta_diff"] <- "Beta_difference"
    names(beta_param)[names(beta_param) == "Delta_beta_diff"] <- "Delta_beta_difference"

    cfile <- file.path(ROOT, "Calibration", "community",
                       paste0(group, "_community_metrics_parameter.csv"))
    if (!file.exists(cfile)) stop("Missing community calibration metrics: ", cfile)
    community <- read.csv(cfile, check.names = FALSE, stringsAsFactors = FALSE)

    all <- merge(beta_param, community, by = param_cols)
    if (!nrow(all)) stop(group, " ", tpl, ": no overlap between beta and community metrics.")

    ww <- weights[weights$TPL == tpl & weights$Group == group, ]
    if (nrow(ww) != 1) stop("Missing/duplicated calibration weight for ", tpl, " ", group)

    all$Beta_component <- all$Beta_difference + all$Delta_beta_difference
    all$Bray_component <- (
      all$TemporalBrayDifference + all$SpatialBrayDifference
    ) / ww$BrayDivisor
    all$Richness_component <- all$RichnessRatioDifference / ww$RichnessDivisor
    all$CalibrationLoss <- all$Beta_component + all$Bray_component + all$Richness_component
    all$BrayDivisor <- ww$BrayDivisor
    all$RichnessDivisor <- ww$RichnessDivisor
    all$EmpiricalBeta_warming <- unname(emp["warming"])
    all$EmpiricalBeta_unwarming <- unname(emp["unwarming"])
    all$Group <- group
    all$TPL <- tpl

    all <- all[order(all$CalibrationLoss, all$Immigration, all$Interaction,
                     all$Growth, all$Noise), ]
    all$Rank <- seq_len(nrow(all))
    top10 <- all[seq_len(min(10L, nrow(all))), ]
    if (nrow(top10) < 10L) stop(group, " ", tpl, ": fewer than 10 candidate parameter sets.")

    outdir <- file.path(ROOT, "Calibration", tpl)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write.csv(all, file.path(outdir, paste0(group, "_all_parameter_loss.csv")),
              row.names = FALSE, quote = FALSE)
    write.csv(top10, file.path(outdir, paste0(group, "_top10_parameters.csv")),
              row.names = FALSE, quote = FALSE)

    message("Selected top 10: ", tpl, " | ", group,
            " | best loss=", signif(top10$CalibrationLoss[1], 6))
  }
}

message("All top-10 sGLV calibration ensembles selected.")
