# ==============================================================================
# NWS SUBGROUP WORKFLOW — STEP 2
# Fig. 2d statistics: treatment-vs-Ctrl Cohen's d + categorical hierarchical LMM
#
# Inputs:
#   ./Taylor_revision_UPARSE_subgroup/
#      tTTPL_phylum_trophic_plot_level_full_precision.rds
#   ./Taylor_revision_pTTPL_subgroup/
#      pTTPL_phylum_trophic_plot_level_full_precision.rds
#
# Main logic
#   A. The first-stage fit scripts already apply the TPL reliability thresholds:
#        Phylum:  n_taxa_used >= 15
#        Trophic: n_taxa_used >= 6
#
#   B. This script applies a SEPARATE Fig. 2d display/representativeness filter:
#        - valid beta in >= 75% of the 48 subplot trajectories
#        - >= 2 valid beta estimates in EVERY treatment cell
#
#   C. For each eligible subgroup and each TPL formulation:
#        - effect size = pooled-SD Cohen's d for each treatment vs Ctrl
#        - formal inference = categorical hierarchical LMM
#             Beta_raw ~ WarmF * PrecipF * ClipF + (1 | block/WholePlotID)
#        - 11 treatment-vs-Ctrl contrasts
#        - Holm adjustment across those 11 contrasts WITHIN
#          TPL_type x Group x Level_type x Subgroup
#
# IMPORTANT
#   - Subgroups are NOT selected because they are significant or have large d.
#   - The Fig. 2d filter is based only on coverage / estimability.
#   - Fig. 2a uses a different estimand (four overall factor effects) and should
#     not be mixed with this treatment-vs-Ctrl analysis.
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

TTTPL_FILE <- "./results/07_tTTPL_subgroup/tTTPL_phylum_trophic_plot_level_full_precision.rds"
PTTPL_FILE <- "./results/08_pTTPL_subgroup/pTTPL_phylum_trophic_plot_level_full_precision.rds"

OUTPUT_ROOT <- "./results/09_subgroup_statistics"
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

TREATMENT_LEVELS <- c(
  "Ctrl", "W", "H", "D", "C",
  "WH", "WD", "WC", "HC", "DC", "WHC", "WDC"
)

TARGET_TREATMENTS <- setdiff(TREATMENT_LEVELS, "Ctrl")

# Display / representativeness filter.
MIN_PROP_VALID_PLOTS <- 0.75
MIN_PER_TREATMENT <- 2L

# "main" reproduces the intended main-figure scope:
#   Archaea/Bacteria/Fungi/Protists -> Phylum
#   Nematoda                         -> Trophic
#
# "all" keeps every fitted phylum/trophic category, including Protist trophic.
FIG2D_SCOPE <- "main"

# Optional plotting. Statistical CSV outputs are created regardless.
MAKE_HEATMAPS <- TRUE

# If TRUE, a paired tTTPL/pTTPL heatmap dataset is restricted to subgroups that
# satisfy the coverage filter for BOTH formulations.
REQUIRE_BOTH_TPL_FOR_PAIRED_HEATMAP <- TRUE

# ------------------------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------------------------

needed <- c("lme4", "lmerTest", "emmeans", "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    "Please install required packages first: ",
    paste(missing, collapse = ", "),
    "\nExample: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))"
  )
}

library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)

# ------------------------------------------------------------------------------
# 2. HELPERS
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

cohen_pooled <- function(x1, x0) {
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]

  n1 <- length(x1)
  n0 <- length(x0)

  if (n1 < 2 || n0 < 2) {
    return(c(
      N_target = n1, N_reference = n0,
      Mean_target = NA, Mean_reference = NA,
      Raw_difference = NA,
      SD_pooled = NA,
      Cohen_d_pooled = NA,
      Hedges_g = NA
    ))
  }

  sp2 <- (
    (n1 - 1) * var(x1) +
    (n0 - 1) * var(x0)
  ) / (n1 + n0 - 2)

  sp <- sqrt(sp2)
  diff <- mean(x1) - mean(x0)

  d <- if (is.finite(sp) && sp > 0) diff / sp else NA_real_

  df <- n1 + n0 - 2
  J <- if (df > 1) 1 - 3 / (4 * df - 1) else NA_real_
  g <- if (is.finite(d) && is.finite(J)) J * d else NA_real_

  c(
    N_target = n1,
    N_reference = n0,
    Mean_target = mean(x1),
    Mean_reference = mean(x0),
    Raw_difference = diff,
    SD_pooled = sp,
    Cohen_d_pooled = d,
    Hedges_g = g
  )
}

sig_label <- function(p) {
  ifelse(
    is.na(p), "",
    ifelse(
      p < 0.001, "***",
      ifelse(p < 0.01, "**",
             ifelse(p < 0.05, "*", ""))
    )
  )
}

scope_keep <- function(group, level_type) {
  if (FIG2D_SCOPE == "all") return(TRUE)

  if (FIG2D_SCOPE == "main") {
    return(
      (group %in% c("Archaea", "Bacteria", "Fungi", "Protists") &&
         level_type == "Phylum") ||
      (group == "Nematoda" && level_type == "Trophic")
    )
  }

  stop("Unknown FIG2D_SCOPE: ", FIG2D_SCOPE)
}

cell_to_treatment <- function(WarmF, PrecipF, ClipF) {
  out <- rep(NA_character_, length(WarmF))

  out[WarmF == "Unwarmed" & PrecipF == "Ambient" & ClipF == "Unclipped"] <- "Ctrl"
  out[WarmF == "Warmed"   & PrecipF == "Ambient" & ClipF == "Unclipped"] <- "W"
  out[WarmF == "Unwarmed" & PrecipF == "Half"    & ClipF == "Unclipped"] <- "H"
  out[WarmF == "Unwarmed" & PrecipF == "Double"  & ClipF == "Unclipped"] <- "D"
  out[WarmF == "Unwarmed" & PrecipF == "Ambient" & ClipF == "Clipped"]   <- "C"
  out[WarmF == "Warmed"   & PrecipF == "Half"    & ClipF == "Unclipped"] <- "WH"
  out[WarmF == "Warmed"   & PrecipF == "Double"  & ClipF == "Unclipped"] <- "WD"
  out[WarmF == "Warmed"   & PrecipF == "Ambient" & ClipF == "Clipped"]   <- "WC"
  out[WarmF == "Unwarmed" & PrecipF == "Half"    & ClipF == "Clipped"]   <- "HC"
  out[WarmF == "Unwarmed" & PrecipF == "Double"  & ClipF == "Clipped"]   <- "DC"
  out[WarmF == "Warmed"   & PrecipF == "Half"    & ClipF == "Clipped"]   <- "WHC"
  out[WarmF == "Warmed"   & PrecipF == "Double"  & ClipF == "Clipped"]   <- "WDC"

  out
}

# ------------------------------------------------------------------------------
# 3. READ AND STANDARDIZE PLOT-LEVEL BETA RESULTS
# ------------------------------------------------------------------------------

if (!file.exists(TTTPL_FILE)) stop("Missing file: ", TTTPL_FILE)
if (!file.exists(PTTPL_FILE)) stop("Missing file: ", PTTPL_FILE)

tt <- readRDS(TTTPL_FILE)
pp <- readRDS(PTTPL_FILE)

required_cols <- c(
  "TPL_type", "Group", "Level_type", "Subgroup", "Plot",
  "block", "WholePlotID",
  "WarmF", "PrecipF", "ClipF", "Treatment",
  "Beta_raw", "SE_Beta_raw"
)

for (nm in required_cols) {
  if (!nm %in% names(tt)) stop("tTTPL file missing column: ", nm)
  if (!nm %in% names(pp)) stop("pTTPL file missing column: ", nm)
}

dat <- rbind(tt[, required_cols], pp[, required_cols])

dat$WarmF <- factor(dat$WarmF,
                    levels = c("Unwarmed", "Warmed"))
dat$PrecipF <- factor(dat$PrecipF,
                      levels = c("Ambient", "Half", "Double"))
dat$ClipF <- factor(dat$ClipF,
                    levels = c("Unclipped", "Clipped"))
dat$Treatment <- factor(dat$Treatment,
                        levels = TREATMENT_LEVELS)
dat$block <- factor(dat$block)
dat$WholePlotID <- factor(dat$WholePlotID)

# Keep only requested figure scope.
keep_scope <- mapply(scope_keep, dat$Group, dat$Level_type)
dat_scope <- dat[keep_scope, , drop = FALSE]

# ------------------------------------------------------------------------------
# 4. COVERAGE / REPRESENTATIVENESS SUMMARY
# ------------------------------------------------------------------------------

keys <- unique(dat_scope[, c("TPL_type", "Group", "Level_type", "Subgroup")])

coverage_rows <- list()
cell_rows <- list()

for (i in seq_len(nrow(keys))) {

  k <- keys[i, , drop = FALSE]
  d <- dat_scope[
    dat_scope$TPL_type == k$TPL_type &
    dat_scope$Group == k$Group &
    dat_scope$Level_type == k$Level_type &
    dat_scope$Subgroup == k$Subgroup,
    ,
    drop = FALSE
  ]

  d_valid <- d[is.finite(d$Beta_raw), , drop = FALSE]

  # The experiment defines 48 subplot trajectories per microbial group.
  n_valid <- length(unique(d_valid$Plot))
  prop_valid <- n_valid / 48

  counts <- setNames(rep(0L, length(TREATMENT_LEVELS)), TREATMENT_LEVELS)
  tab <- table(as.character(d_valid$Treatment))
  counts[names(tab)] <- as.integer(tab)

  min_per_treatment <- min(counts)
  eligible <- (
    prop_valid >= MIN_PROP_VALID_PLOTS &&
    min_per_treatment >= MIN_PER_TREATMENT
  )

  coverage_rows[[length(coverage_rows) + 1]] <- data.frame(
    TPL_type = k$TPL_type,
    Group = k$Group,
    Level_type = k$Level_type,
    Subgroup = k$Subgroup,
    N_valid_plots = n_valid,
    Proportion_valid = prop_valid,
    Min_n_per_treatment = min_per_treatment,
    Eligible_for_Fig2d = eligible,
    stringsAsFactors = FALSE
  )

  for (tr in TREATMENT_LEVELS) {
    cell_rows[[length(cell_rows) + 1]] <- data.frame(
      TPL_type = k$TPL_type,
      Group = k$Group,
      Level_type = k$Level_type,
      Subgroup = k$Subgroup,
      Treatment = tr,
      N_valid = counts[[tr]],
      stringsAsFactors = FALSE
    )
  }
}

coverage <- do.call(rbind, coverage_rows)
cell_counts <- do.call(rbind, cell_rows)

# Eligibility in BOTH formulations, useful when plotting beta_t and beta_p side-by-side.
pair_key <- unique(coverage[, c("Group", "Level_type", "Subgroup")])
both_rows <- list()

for (i in seq_len(nrow(pair_key))) {
  k <- pair_key[i, , drop = FALSE]
  z <- coverage[
    coverage$Group == k$Group &
    coverage$Level_type == k$Level_type &
    coverage$Subgroup == k$Subgroup,
    ,
    drop = FALSE
  ]

  has_t <- any(z$TPL_type == "tTTPL" & z$Eligible_for_Fig2d)
  has_p <- any(z$TPL_type == "pTTPL" & z$Eligible_for_Fig2d)

  both_rows[[length(both_rows) + 1]] <- data.frame(
    Group = k$Group,
    Level_type = k$Level_type,
    Subgroup = k$Subgroup,
    Eligible_tTTPL = has_t,
    Eligible_pTTPL = has_p,
    Eligible_both_TPL = has_t & has_p,
    stringsAsFactors = FALSE
  )
}

elig_both <- do.call(rbind, both_rows)

write_csv_full_precision(
  cell_counts,
  file.path(OUTPUT_ROOT, "01_subgroup_cell_counts.csv")
)
write_csv_full_precision(
  coverage,
  file.path(OUTPUT_ROOT, "01b_subgroup_eligibility_by_TPL.csv")
)
write_csv_full_precision(
  elig_both,
  file.path(OUTPUT_ROOT, "01c_subgroup_eligibility_both_TPL.csv")
)

# ------------------------------------------------------------------------------
# 5. COHEN'S d + CATEGORICAL HIERARCHICAL LMM
# ------------------------------------------------------------------------------

eligible_keys <- coverage[
  coverage$Eligible_for_Fig2d,
  c("TPL_type", "Group", "Level_type", "Subgroup"),
  drop = FALSE
]

lmm_rows <- list()
d_rows <- list()
model_qc_rows <- list()

for (i in seq_len(nrow(eligible_keys))) {

  k <- eligible_keys[i, , drop = FALSE]

  d <- dat_scope[
    dat_scope$TPL_type == k$TPL_type &
    dat_scope$Group == k$Group &
    dat_scope$Level_type == k$Level_type &
    dat_scope$Subgroup == k$Subgroup &
    is.finite(dat_scope$Beta_raw),
    ,
    drop = FALSE
  ]

  # Reset factor levels after subset.
  d$WarmF <- factor(d$WarmF, levels = c("Unwarmed", "Warmed"))
  d$PrecipF <- factor(d$PrecipF, levels = c("Ambient", "Half", "Double"))
  d$ClipF <- factor(d$ClipF, levels = c("Unclipped", "Clipped"))
  d$Treatment <- factor(d$Treatment, levels = TREATMENT_LEVELS)
  d$block <- factor(d$block)
  d$WholePlotID <- factor(d$WholePlotID)

  fit_error <- NA_character_
  fit <- tryCatch(
    lmer(
      Beta_raw ~ WarmF * PrecipF * ClipF +
        (1 | block/WholePlotID),
      data = d,
      REML = TRUE
    ),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    model_qc_rows[[length(model_qc_rows) + 1]] <- data.frame(
      TPL_type = k$TPL_type,
      Group = k$Group,
      Level_type = k$Level_type,
      Subgroup = k$Subgroup,
      N = nrow(d),
      Singular = NA,
      Model_error = fit_error,
      stringsAsFactors = FALSE
    )
    next
  }

  singular_flag <- isSingular(fit, tol = 1e-4)

  model_qc_rows[[length(model_qc_rows) + 1]] <- data.frame(
    TPL_type = k$TPL_type,
    Group = k$Group,
    Level_type = k$Level_type,
    Subgroup = k$Subgroup,
    N = nrow(d),
    Singular = singular_flag,
    Model_error = NA_character_,
    stringsAsFactors = FALSE
  )

  emm <- emmeans(
    fit,
    ~ WarmF * PrecipF * ClipF
  )

  emm_df <- as.data.frame(emm)
  emm_df$Treatment <- cell_to_treatment(
    as.character(emm_df$WarmF),
    as.character(emm_df$PrecipF),
    as.character(emm_df$ClipF)
  )

  if (any(is.na(emm_df$Treatment))) {
    stop("Could not map at least one emmeans cell to a treatment label for ",
         paste(unlist(k), collapse = " / "))
  }

  contrast_list <- list()

  for (tr in TARGET_TREATMENTS) {
    v <- rep(0, nrow(emm_df))
    v[emm_df$Treatment == tr] <- 1
    v[emm_df$Treatment == "Ctrl"] <- -1
    contrast_list[[paste0(tr, "-Ctrl")]] <- v
  }

  cc <- contrast(emm, method = contrast_list, adjust = "none")
  cs <- as.data.frame(summary(cc, infer = c(TRUE, TRUE)))

  # Rename CI columns robustly.
  low_col <- intersect(c("lower.CL", "asymp.LCL"), names(cs))
  high_col <- intersect(c("upper.CL", "asymp.UCL"), names(cs))
  t_col <- intersect(c("t.ratio", "z.ratio"), names(cs))

  if (length(low_col) != 1 || length(high_col) != 1 || length(t_col) != 1) {
    stop("Unexpected emmeans summary columns.")
  }

  cs$Treatment <- sub("-Ctrl$", "", cs$contrast)
  cs$P_Holm_11_vsCtrl <- p.adjust(cs$p.value, method = "holm")
  cs$Significance_Holm <- sig_label(cs$P_Holm_11_vsCtrl)

  for (j in seq_len(nrow(cs))) {
    lmm_rows[[length(lmm_rows) + 1]] <- data.frame(
      TPL_type = k$TPL_type,
      Group = k$Group,
      Level_type = k$Level_type,
      Subgroup = k$Subgroup,
      Treatment = cs$Treatment[j],
      Reference = "Ctrl",
      Contrast = cs$contrast[j],
      Estimate = cs$estimate[j],
      SE = cs$SE[j],
      df = if ("df" %in% names(cs)) cs$df[j] else NA_real_,
      CI_low = cs[[low_col]][j],
      CI_high = cs[[high_col]][j],
      Test_statistic = cs[[t_col]][j],
      P_raw = cs$p.value[j],
      P_Holm_11_vsCtrl = cs$P_Holm_11_vsCtrl[j],
      Significance_Holm = cs$Significance_Holm[j],
      Significant_Holm_11_vsCtrl =
        is.finite(cs$P_Holm_11_vsCtrl[j]) &&
        cs$P_Holm_11_vsCtrl[j] < 0.05,
      Model_singular = singular_flag,
      stringsAsFactors = FALSE
    )
  }

  # Descriptive standardized effect sizes from raw subplot beta values.
  ctrl <- d$Beta_raw[d$Treatment == "Ctrl"]

  for (tr in TARGET_TREATMENTS) {
    target <- d$Beta_raw[d$Treatment == tr]
    es <- cohen_pooled(target, ctrl)

    d_rows[[length(d_rows) + 1]] <- data.frame(
      TPL_type = k$TPL_type,
      Group = k$Group,
      Level_type = k$Level_type,
      Subgroup = k$Subgroup,
      Treatment = tr,
      Reference = "Ctrl",
      N_target = unname(es["N_target"]),
      N_reference = unname(es["N_reference"]),
      Mean_target = unname(es["Mean_target"]),
      Mean_reference = unname(es["Mean_reference"]),
      Raw_difference = unname(es["Raw_difference"]),
      SD_pooled = unname(es["SD_pooled"]),
      Cohen_d_pooled = unname(es["Cohen_d_pooled"]),
      Hedges_g = unname(es["Hedges_g"]),
      stringsAsFactors = FALSE
    )
  }
}

lmm <- if (length(lmm_rows) > 0) do.call(rbind, lmm_rows) else data.frame()
eff <- if (length(d_rows) > 0) do.call(rbind, d_rows) else data.frame()
model_qc <- if (length(model_qc_rows) > 0) do.call(rbind, model_qc_rows) else data.frame()

if (nrow(lmm) == 0 || nrow(eff) == 0) {
  stop(
    "No Fig. 2d subgroup statistics were generated. ",
    "Inspect 01b_subgroup_eligibility_by_TPL.csv first."
  )
}

write_csv_full_precision(
  lmm,
  file.path(OUTPUT_ROOT, "02_subgroup_11_vsCtrl_LMM_Holm.csv")
)
write_csv_full_precision(
  eff,
  file.path(OUTPUT_ROOT, "03_subgroup_11_vsCtrl_Cohen_d_Hedges_g.csv")
)
write_csv_full_precision(
  model_qc,
  file.path(OUTPUT_ROOT, "03b_subgroup_LMM_QC.csv")
)

# ------------------------------------------------------------------------------
# 6. MERGED HEATMAP-READY TABLE
# ------------------------------------------------------------------------------

merge_keys <- c(
  "TPL_type", "Group", "Level_type",
  "Subgroup", "Treatment", "Reference"
)

heat <- merge(
  eff, lmm,
  by = merge_keys,
  all = TRUE,
  sort = FALSE
)

heat$Sig <- sig_label(heat$P_Holm_11_vsCtrl)

# Preserve the intended treatment order.
heat$Treatment <- factor(
  heat$Treatment,
  levels = TARGET_TREATMENTS
)

# Attach display eligibility information.
heat <- merge(
  heat,
  coverage[, c(
    "TPL_type", "Group", "Level_type", "Subgroup",
    "N_valid_plots", "Proportion_valid",
    "Min_n_per_treatment", "Eligible_for_Fig2d"
  )],
  by = c("TPL_type", "Group", "Level_type", "Subgroup"),
  all.x = TRUE,
  sort = FALSE
)

heat <- merge(
  heat,
  elig_both,
  by = c("Group", "Level_type", "Subgroup"),
  all.x = TRUE,
  sort = FALSE
)

heat$Treatment <- factor(
  heat$Treatment,
  levels = TARGET_TREATMENTS
)

write_csv_full_precision(
  heat,
  file.path(OUTPUT_ROOT, "04_Fig2d_heatmap_ready_major_by_TPL.csv")
)
saveRDS(
  heat,
  file.path(OUTPUT_ROOT, "04_Fig2d_heatmap_ready_major_by_TPL.rds")
)

if (REQUIRE_BOTH_TPL_FOR_PAIRED_HEATMAP) {
  heat_both <- heat[heat$Eligible_both_TPL %in% TRUE, , drop = FALSE]
} else {
  heat_both <- heat
}

write_csv_full_precision(
  heat_both,
  file.path(OUTPUT_ROOT, "05_Fig2d_heatmap_ready_major_both_TPL.csv")
)

# ------------------------------------------------------------------------------
# 7. OPTIONAL HEATMAPS
# ------------------------------------------------------------------------------

if (MAKE_HEATMAPS) {

  plot_one <- function(z, tpl_type, filename) {

    z <- z[z$TPL_type == tpl_type, , drop = FALSE]
    if (nrow(z) == 0) return(invisible(NULL))

    # Order subgroups within each microbial group by median absolute d,
    # ONLY for visual ordering, not for inclusion/exclusion.
    ord <- aggregate(
      abs(Cohen_d_pooled) ~ Group + Subgroup,
      data = z,
      FUN = function(x) median(x, na.rm = TRUE)
    )
    ord <- ord[order(ord$Group, -ord[[3]]), ]

    # Create unique y labels so identical subgroup names in different groups
    # cannot collide internally; facet strip still gives the microbial group.
    z$Subgroup_plot <- paste(z$Group, z$Subgroup, sep = "___")

    level_order <- unlist(lapply(unique(ord$Group), function(g) {
      paste(
        g,
        ord$Subgroup[ord$Group == g],
        sep = "___"
      )
    }))

    z$Subgroup_plot <- factor(
      z$Subgroup_plot,
      levels = rev(unique(level_order))
    )

    label_map <- setNames(
      sub("^.*___", "", levels(z$Subgroup_plot)),
      levels(z$Subgroup_plot)
    )

    p <- ggplot(
      z,
      aes(
        x = Treatment,
        y = Subgroup_plot,
        fill = Cohen_d_pooled
      )
    ) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(
        aes(label = Sig),
        size = 3.2,
        fontface = "bold"
      ) +
      scale_y_discrete(labels = label_map) +
      scale_fill_gradient2(
        midpoint = 0,
        name = "Cohen's d"
      ) +
      facet_grid(
        Group ~ .,
        scales = "free_y",
        space = "free_y"
      ) +
      labs(
        x = NULL,
        y = NULL,
        title = tpl_type
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
          angle = 45, hjust = 1, vjust = 1,
          color = "black"
        ),
        axis.text.y = element_text(color = "black"),
        strip.text.y = element_text(
          angle = 0, face = "bold"
        ),
        strip.background = element_blank()
      )

    ggsave(
      file.path(OUTPUT_ROOT, filename),
      p,
      width = 8,
      height = max(4, 0.25 * length(unique(z$Subgroup_plot)) + 1.5)
    )

    invisible(p)
  }

  plot_one(
    heat,
    "tTTPL",
    "Fig2d_tTTPL_major_heatmap.pdf"
  )

  plot_one(
    heat,
    "pTTPL",
    "FigS_pTTPL_major_heatmap.pdf"
  )
}

# ------------------------------------------------------------------------------
# 8. CONSOLE SUMMARY
# ------------------------------------------------------------------------------

message("\nDisplay filter:")
message("  Proportion valid plots >= ", MIN_PROP_VALID_PLOTS)
message("  Minimum valid beta per treatment >= ", MIN_PER_TREATMENT)
message("  Scope = ", FIG2D_SCOPE)

message("\nEligible subgroup counts by TPL type:")
print(
  aggregate(
    Eligible_for_Fig2d ~ TPL_type,
    data = coverage,
    FUN = sum
  )
)

message("\nSaved Fig. 2d statistics to: ", OUTPUT_ROOT)
message("Done.")
