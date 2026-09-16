# ==============================================================================
# pTTPL tip-shuffle null analysis — end-to-end
# ==============================================================================
# This script replaces three legacy files:
#   Calculate_null_MPDik.R
#   Taylor_law_power_pTTPL_test_null_model.R
#   2.Taylor_law_power_compare_null.R
#
# For each soil-biotic group it:
#   1) tip-shuffles the empirical phylogeny 100 times by default;
#   2) recalculates MPD_i,k while retaining empirical abundance trajectories;
#   3) refits beta_p separately for every subplot under every shuffle;
#   4) averages the 100 null beta_p values within each subplot;
#   5) compares observed and matched null expectations with a two-sided paired
#      Wilcoxon signed-rank test;
#   6) Holm-adjusts the five group-level P values.
#
# The null realizations are Monte Carlo phylogenetic randomizations and are NOT
# treated as independent biological replicates.
# ==============================================================================

rm(list = ls())

PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
setwd(PROJECT_DIR)

PHYLO_DATA_DIR <- Sys.getenv(
  "TPL_PHYLO_DATA_DIR",
  unset = file.path(PROJECT_DIR, "data")
)
META_FILE <- Sys.getenv(
  "TPL_METADATA_FILE",
  unset = file.path(PROJECT_DIR, "data", "prepared", "environment_metadata.csv")
)
OBSERVED_FILE <- file.path(
  PROJECT_DIR,
  "results", "04_pTTPL",
  "pTTPL_all_groups_plot_level_full_precision.csv"
)
OUT_DIR <- file.path(PROJECT_DIR, "results", "06_tipshuffle_null")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_SHUFFLES <- as.integer(Sys.getenv("TPL_TIPSHUFFLE_N", unset = "100"))
if (!is.finite(N_SHUFFLES) || N_SHUFFLES < 1L) stop("TPL_TIPSHUFFLE_N must be >= 1.")

# The legacy analysis did not record a fixed random seed. To preserve that
# behavior, no seed is imposed unless TPL_TIPSHUFFLE_SEED is explicitly set.
seed_text <- Sys.getenv("TPL_TIPSHUFFLE_SEED", unset = "")
if (nzchar(seed_text)) set.seed(as.integer(seed_text))

group_text <- Sys.getenv(
  "TPL_TIPSHUFFLE_GROUPS",
  unset = "Archaea,Bacteria,Fungi,Protists,Nematoda"
)
GROUPS <- trimws(strsplit(group_text, ",", fixed = TRUE)[[1]])

required_packages <- c("iCAMP", "ape", "bigmemory")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R package(s): ", paste(missing_packages, collapse = ", "))
}

read_numeric_matrix <- function(file) {
  if (!file.exists(file)) stop("Input file not found: ", file)
  x <- read.csv(file, header = TRUE, row.names = 1, check.names = FALSE,
                stringsAsFactors = FALSE)
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) stop("Non-numeric/NA values in: ", file)
  x
}

find_input <- function(...) {
  candidates <- c(...)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("None of the candidate input files exists:\n", paste(candidates, collapse = "\n"))
  }
  hit[1]
}

mpdi_big <- function(comm, pdist.out, sp.limit = 10000L) {
  # comm: samples x taxa
  ids <- match(colnames(comm), pdist.out$tip.label)
  if (anyNA(ids)) stop("Some taxa could not be matched to the shuffled distance matrix.")

  rs <- rowSums(comm)
  if (any(!is.finite(rs)) || any(rs <= 0)) stop("Invalid sample totals in community matrix.")
  comm <- comm / rs

  pd <- try(
    bigmemory::attach.big.matrix(
      dget(file.path(pdist.out$pd.wd, pdist.out$pd.file))
    ),
    silent = TRUE
  )
  if (inherits(pd, "try-error")) {
    pd <- bigmemory::attach.big.matrix(file.path(pdist.out$pd.wd, pdist.out$pd.file))
  }

  if (nrow(pd) <= sp.limit) {
    comd <- comm %*% pd[ids, ids]
  } else {
    nsp <- ncol(comm)
    n_int <- max(1L, floor(sp.limit * sp.limit / nsp))
    starts <- seq.int(1L, nsp, by = n_int)
    comd <- matrix(0, nrow = nrow(comm), ncol = nsp)
    for (from in starts) {
      to <- min(from + n_int - 1L, nsp)
      comd[, from:to] <- comm %*% pd[ids, ids[from:to], drop = FALSE]
    }
  }

  out <- comd * comm
  rownames(out) <- rownames(comm)
  colnames(out) <- colnames(comm)
  out
}

fit_beta <- function(taxon_by_year) {
  mu <- rowMeans(taxon_by_year)
  vv <- apply(taxon_by_year, 1L, var)
  keep <- is.finite(mu) & is.finite(vv) & mu > 0 & vv > 0
  if (sum(keep) < 3L) return(c(Intercept = NA_real_, Beta = NA_real_, R2 = NA_real_, P = NA_real_))
  dat <- data.frame(LogMean = log10(mu[keep]), LogVar = log10(vv[keep]))
  fit <- lm(LogVar ~ LogMean, data = dat)
  sm <- summary(fit)
  c(
    Intercept = unname(coef(fit)[1]),
    Beta = unname(coef(fit)[2]),
    R2 = unname(sm$r.squared),
    P = unname(sm$coefficients[2, 4])
  )
}

group_config <- list(
  Archaea = list(
    comm = "resample7842_otutab_Archaea_2010to2020.csv",
    tree = "rooted_tree_otus_Archaea_2010to2020_resample7842.nwk"
  ),
  Bacteria = list(
    comm = "resample26890_otutab_Bacteria_2010to2020.csv",
    tree = "rooted_tree_otus_Bacteria_2010to2020_resample26890.nwk"
  ),
  Fungi = list(
    comm = "resample541_otutab_Fungi_2010to2020.csv",
    tree = "rooted_tree_otus_Fungi_2010to2020_resample541.nwk"
  ),
  Protists = list(
    comm = "resample1064_otutab_Protists_2010to2020.csv",
    tree = "rooted_tree_otus_Protists_2010to2020_resample1064.nwk"
  ),
  Nematoda = list(
    comm = "resample296_otutab_Nematoda_2010to2020.csv",
    tree = "rooted_tree_otus_Nematoda_2010to2020_resample296.nwk"
  )
)

if (!all(GROUPS %in% names(group_config))) {
  stop("Unknown TPL_TIPSHUFFLE_GROUPS value(s): ",
       paste(setdiff(GROUPS, names(group_config)), collapse = ", "))
}
if (!file.exists(META_FILE)) stop("Metadata file not found: ", META_FILE)
if (!file.exists(OBSERVED_FILE)) {
  stop("Observed full-precision pTTPL file not found: ", OBSERVED_FILE,
       "\nRun 04_pTTPL_fit_full_precision.R first.")
}

meta <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)
needed_meta <- c("Sample_name", "plot", "year")
miss <- setdiff(needed_meta, names(meta))
if (length(miss)) stop("Metadata missing column(s): ", paste(miss, collapse = ", "))
meta$Sample_name <- as.character(meta$Sample_name)
meta$plot <- as.character(meta$plot)
meta$year <- as.numeric(meta$year)
meta <- meta[meta$year %in% 2010:2020, , drop = FALSE]

null_rows <- list()
row_id <- 1L

for (group in GROUPS) {
  message("\n=== Tip-shuffle null: ", group, " ===")
  cfg <- group_config[[group]]

  comm_file <- find_input(
    file.path(PHYLO_DATA_DIR, cfg$comm),
    file.path(PHYLO_DATA_DIR, cfg$comm)
  )
  tree_file <- find_input(
    file.path(PHYLO_DATA_DIR, cfg$tree),
    file.path(PHYLO_DATA_DIR, cfg$tree)
  )

  comm <- read_numeric_matrix(comm_file)  # taxa x samples
  tree <- ape::read.tree(tree_file)

  common_taxa <- intersect(rownames(comm), tree$tip.label)
  if (length(common_taxa) < 3L) stop(group, ": fewer than 3 shared taxa.")
  comm <- comm[rownames(comm) %in% common_taxa, , drop = FALSE]
  tree <- ape::keep.tip(tree, common_taxa)

  md_group <- meta[meta$Sample_name %in% colnames(comm), , drop = FALSE]
  plot_list <- unique(md_group$plot)
  if (!length(plot_list)) stop(group, ": no sample names match metadata.")
  comm <- comm[, unique(md_group$Sample_name), drop = FALSE]
  comm_samples <- t(comm)

  for (shuffle_id in seq_len(N_SHUFFLES)) {
    message(group, ": shuffle ", shuffle_id, "/", N_SHUFFLES)

    shuffled_tree <- tree
    shuffled_tree$tip.label <- sample(tree$tip.label, replace = FALSE)

    tmp <- tempfile(pattern = paste0("tipshuffle_", group, "_", shuffle_id, "_"))
    dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
    pdist.out <- iCAMP::pdist.big(shuffled_tree, wd = tmp)
    community_mpd <- iCAMP::mpd.big(
      comm_samples,
      pd.desc = pdist.out$pd.file,
      pd.spname = pdist.out$tip.label,
      pd.wd = pdist.out$pd.wd,
      spname.check = TRUE
    )
    each_taxon <- mpdi_big(comm_samples, pdist.out)
    mpdik <- sweep(each_taxon, 1L, community_mpd, "/")
    if (any(!is.finite(mpdik))) stop(group, ": non-finite MPD_i,k in shuffle ", shuffle_id)

    for (plot_id in plot_list) {
      md <- md_group[md_group$plot == plot_id, , drop = FALSE]
      md <- md[order(md$year), , drop = FALSE]
      ids <- md$Sample_name
      if (!all(ids %in% rownames(mpdik))) {
        stop(group, " / ", plot_id, ": unmatched annual samples during null analysis.")
      }
      fit <- fit_beta(t(mpdik[ids, , drop = FALSE]))
      null_rows[[row_id]] <- data.frame(
        Group = group,
        Null = shuffle_id,
        Plot = plot_id,
        Intercept = fit[["Intercept"]],
        Beta = fit[["Beta"]],
        R2 = fit[["R2"]],
        P = fit[["P"]],
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }

    unlink(tmp, recursive = TRUE, force = TRUE)
  }
}

null_beta <- do.call(rbind, null_rows)
write.csv(null_beta, file.path(OUT_DIR, "tipshuffle_null_beta_all.csv"),
          row.names = FALSE, quote = FALSE)

observed <- read.csv(OBSERVED_FILE, check.names = FALSE, stringsAsFactors = FALSE)
needed_obs <- c("Group", "Plot", "Beta_raw")
miss <- setdiff(needed_obs, names(observed))
if (length(miss)) stop("Observed pTTPL table missing: ", paste(miss, collapse = ", "))
observed <- observed[observed$Group %in% GROUPS, needed_obs, drop = FALSE]
names(observed)[names(observed) == "Beta_raw"] <- "Beta_observed"

null_mean <- aggregate(Beta ~ Group + Plot, null_beta, mean, na.rm = TRUE)
names(null_mean)[names(null_mean) == "Beta"] <- "Beta_null_mean"
paired <- merge(observed, null_mean, by = c("Group", "Plot"), all = FALSE)
paired$Delta_observed_minus_null <- paired$Beta_observed - paired$Beta_null_mean
write.csv(paired, file.path(OUT_DIR, "tipshuffle_observed_vs_null_by_plot.csv"),
          row.names = FALSE, quote = FALSE)

result_rows <- lapply(GROUPS, function(group) {
  d <- paired[paired$Group == group, , drop = FALSE]
  if (nrow(d) < 2L) stop(group, ": fewer than two matched plot pairs.")
  wt <- wilcox.test(
    d$Beta_observed,
    d$Beta_null_mean,
    paired = TRUE,
    alternative = "two.sided",
    exact = FALSE
  )
  data.frame(
    Group = group,
    N_matched_plots = nrow(d),
    Mean_observed = mean(d$Beta_observed, na.rm = TRUE),
    Mean_null_expectation = mean(d$Beta_null_mean, na.rm = TRUE),
    Mean_difference = mean(d$Delta_observed_minus_null, na.rm = TRUE),
    Median_difference = median(d$Delta_observed_minus_null, na.rm = TRUE),
    Wilcoxon_V = unname(wt$statistic),
    P_raw = wt$p.value,
    stringsAsFactors = FALSE
  )
})
paired_tests <- do.call(rbind, result_rows)
paired_tests$P_Holm_5groups <- p.adjust(paired_tests$P_raw, method = "holm")
write.csv(paired_tests, file.path(OUT_DIR, "tipshuffle_paired_tests.csv"),
          row.names = FALSE, quote = FALSE)

message("\nTip-shuffle analysis completed. Outputs: ", OUT_DIR)
print(paired_tests)
