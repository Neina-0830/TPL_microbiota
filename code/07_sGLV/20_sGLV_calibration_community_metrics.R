# ==============================================================================
# 20_sGLV_calibration_community_metrics.R
# Community-level calibration discrepancies for the sGLV model.
#
# External empirical inputs are read directly from ./data:
#   - resampled taxon tables for the five soil-biotic groups
#   - summary of environmental data2010-2020.csv (sample/treatment metadata)
#
# The script reconstructs the empirical calibration targets used by the supplied
# legacy workflow rather than requiring additional precomputed target files:
#   1) temporal Bray-Curtis dissimilarity, averaged across plots within warming
#      and unwarming conditions for each pair of years;
#   2) warming-vs-unwarming spatial Bray-Curtis dissimilarity, averaged across
#      matched precipitation x clipping backgrounds and blocks for each year;
#   3) warming/unwarming richness ratio, averaged across the same matched cells
#      for each year.
#
# Simulated calibration trajectories are read from ./results/17_sGLV/Resample.
# ============================================================================== 

rm(list = ls())

PROJECT_DIR <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()),
                             winslash = "/", mustWork = TRUE)
DATA_DIR <- file.path(PROJECT_DIR, "data")
ROOT <- normalizePath(Sys.getenv("TPL_SGLV_WORK_DIR",
                                 unset = file.path(PROJECT_DIR, "results", "17_sGLV")),
                      winslash = "/", mustWork = FALSE)
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

if (!requireNamespace("vegan", quietly = TRUE)) stop("Install R package 'vegan'.")

CONFIG_FILE <- file.path(PROJECT_DIR, "code", "07_sGLV", "sGLV_group_config.csv")
cfg <- read.csv(CONFIG_FILE, check.names = FALSE, stringsAsFactors = FALSE)

groups_env <- Sys.getenv("TPL_SGLV_GROUPS", unset = "all")
GROUPS <- if (tolower(groups_env) == "all") cfg$Group else trimws(strsplit(groups_env, ",")[[1]])

RAW_COMM_FILES <- c(
  Archaea  = "resample7842_otutab_Archaea_2010to2020.csv",
  Bacteria = "resample26890_otutab_Bacteria_2010to2020.csv",
  Fungi_18S = "resample541_otutab_Fungi_2010to2020.csv",
  Protists = "resample1064_otutab_Protists_2010to2020.csv",
  Nematoda = "resample296_otutab_Nematoda_2010to2020.csv"
)
META_FILE <- file.path(DATA_DIR, "summary of environmental data2010-2020.csv")

parse_calibration_filename <- function(file_name) {
  pat <- paste0(
    "^Time_series_from_SLM_",
    "immigration_d([^_]+)_interaction_sd([^_]+)_growth_g([^_]+)_",
    "seed([^_]+)_noise([^\\.]+)\\.csv$"
  )
  hit <- regmatches(basename(file_name), regexec(pat, basename(file_name)))[[1]]
  if (length(hit) != 6) stop("Cannot parse calibration filename: ", file_name)
  data.frame(
    Immigration = as.numeric(hit[2]),
    Interaction = as.numeric(hit[3]),
    Growth = as.numeric(hit[4]),
    Seed = as.character(hit[5]),
    Noise = as.numeric(hit[6]),
    stringsAsFactors = FALSE
  )
}

read_counts <- function(file) {
  x <- read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 11) stop("Fewer than 11 annual columns in: ", file)
  x <- as.matrix(x[, seq_len(11), drop = FALSE])
  storage.mode(x) <- "double"
  x
}

# Canonical sample metadata from the uploaded environmental table.
read_sample_metadata <- function() {
  if (!file.exists(META_FILE)) stop("Missing empirical metadata: ", META_FILE)
  md <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)
  names(md)[1] <- "Sample_name"
  required <- c("Sample_name", "plot", "block", "year", "Warm", "Precipitation", "Clip")
  miss <- setdiff(required, names(md))
  if (length(miss)) stop("Metadata is missing: ", paste(miss, collapse = ", "))
  md$Sample_name <- as.character(md$Sample_name)
  md$plot <- as.character(md$plot)
  md$block <- as.character(md$block)
  md$year <- as.integer(md$year)
  md$condition <- ifelse(tolower(as.character(md$Warm)) == "warming" |
                           suppressWarnings(as.numeric(md$warm)) == 1,
                         "warming", "unwarming")
  # Historical 'background' matched warming/unwarming within the same
  # precipitation x clipping context. Reconstruct that pairing explicitly.
  md$background <- interaction(md$Precipitation, md$Clip, drop = TRUE, sep = "__")
  md
}

read_empirical_counts <- function(group, md) {
  f <- file.path(DATA_DIR, unname(RAW_COMM_FILES[group]))
  if (!file.exists(f)) stop("Missing empirical community table: ", f)
  x <- read.csv(f, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  missing_samples <- setdiff(md$Sample_name, colnames(x))
  if (length(missing_samples)) {
    stop(group, ": empirical table is missing sample(s), e.g. ",
         paste(head(missing_samples, 10), collapse = ", "))
  }
  x[, md$Sample_name, drop = FALSE]
}

pair_key <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  ifelse(a <= b, paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}

bray_from_sample_matrix <- function(sample_by_taxa, labels) {
  d <- as.matrix(vegan::vegdist(sample_by_taxa, method = "bray"))
  ind <- which(upper.tri(d), arr.ind = TRUE)
  data.frame(
    name1 = labels[ind[, 1]],
    name2 = labels[ind[, 2]],
    bray = d[ind],
    stringsAsFactors = FALSE
  )
}

# Match the legacy empirical target: calculate temporal beta diversity within
# each plot, convert sample labels to years, then average a given year-pair
# across all plots in the same warming condition.
empirical_temporal_target <- function(counts, md, condition) {
  sub <- md[md$condition == condition, , drop = FALSE]
  plots <- unique(sub$plot)
  tabs <- lapply(plots, function(p) {
    z <- sub[sub$plot == p, , drop = FALSE]
    z <- z[order(z$year), , drop = FALSE]
    if (!identical(z$year, 2010:2020)) stop(condition, " / ", p, ": incomplete 2010-2020 series")
    mat <- t(counts[, z$Sample_name, drop = FALSE])
    bray_from_sample_matrix(mat, paste0("Y", z$year))
  })
  all <- do.call(rbind, tabs)
  all$key <- pair_key(all$name1, all$name2)
  agg <- aggregate(bray ~ key, all, mean, na.rm = TRUE)
  split_key <- strsplit(agg$key, "\\|\\|")
  agg$name1 <- vapply(split_key, `[`, character(1), 1)
  agg$name2 <- vapply(split_key, `[`, character(1), 2)
  agg[, c("name1", "name2", "bray")]
}

# Match the legacy spatial target: within each precipitation x clipping
# background, block and year, compare the corresponding warmed and unwarmed
# sample; then average that distance across backgrounds/blocks for each year.
empirical_spatial_target <- function(counts, md) {
  years <- sort(unique(md$year))
  rows <- list(); k <- 0L
  for (bg in unique(md$background)) {
    for (bl in unique(md$block)) {
      for (yr in years) {
        z <- md[md$background == bg & md$block == bl & md$year == yr, , drop = FALSE]
        w <- z$Sample_name[z$condition == "warming"]
        u <- z$Sample_name[z$condition == "unwarming"]
        if (length(w) != 1L || length(u) != 1L) {
          stop("Could not form a unique warming/unwarming pair for background=", bg,
               ", block=", bl, ", year=", yr)
        }
        m <- rbind(
          warming = counts[, w],
          unwarming = counts[, u]
        )
        d <- as.matrix(vegan::vegdist(m, method = "bray"))[1, 2]
        k <- k + 1L
        rows[[k]] <- data.frame(
          name1 = paste0("Y", yr, "_warming"),
          name2 = paste0("Y", yr, "_unwarming"),
          bray = d,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  all <- do.call(rbind, rows)
  aggregate(bray ~ name1 + name2, all, mean, na.rm = TRUE)
}

# Match the legacy richness target: calculate warmed/unwarmed richness ratio in
# every matched background x block x year cell and average across cells by year.
empirical_richness_ratio_target <- function(counts, md) {
  richness <- colSums(counts > 0)
  years <- sort(unique(md$year))
  out <- numeric(length(years)); names(out) <- paste0("Y", years)
  for (ii in seq_along(years)) {
    yr <- years[ii]
    vals <- numeric()
    for (bg in unique(md$background)) {
      for (bl in unique(md$block)) {
        z <- md[md$background == bg & md$block == bl & md$year == yr, , drop = FALSE]
        w <- z$Sample_name[z$condition == "warming"]
        u <- z$Sample_name[z$condition == "unwarming"]
        if (length(w) != 1L || length(u) != 1L) next
        vals <- c(vals, richness[w] / richness[u])
      }
    }
    if (!length(vals)) stop("No empirical richness ratios for year ", yr)
    out[ii] <- mean(vals, na.rm = TRUE)
  }
  out
}

mse_pair_tables <- function(emp, sim, label) {
  emp$key <- pair_key(emp$name1, emp$name2)
  sim$key <- pair_key(sim$name1, sim$name2)
  m <- merge(emp[, c("key", "bray")], sim[, c("key", "bray")],
             by = "key", suffixes = c("_emp", "_sim"))
  if (!nrow(m)) stop(label, ": no matching pair labels.")
  mean((m$bray_emp - m$bray_sim)^2, na.rm = TRUE)
}

bray_pairs <- function(counts_taxa_by_year) {
  years <- colnames(counts_taxa_by_year)
  mat <- t(counts_taxa_by_year)
  bray_from_sample_matrix(mat, years)
}

spatial_pairs <- function(warm_counts, unwarm_counts) {
  years <- intersect(colnames(warm_counts), colnames(unwarm_counts))
  years <- years[seq_len(min(11L, length(years)))]
  taxa <- union(rownames(warm_counts), rownames(unwarm_counts))
  w <- matrix(0, nrow = length(taxa), ncol = length(years), dimnames = list(taxa, years))
  u <- w
  w[rownames(warm_counts), ] <- warm_counts[, years, drop = FALSE]
  u[rownames(unwarm_counts), ] <- unwarm_counts[, years, drop = FALSE]
  do.call(rbind, lapply(years, function(y) {
    m <- rbind(warming = w[, y], unwarming = u[, y])
    d <- as.matrix(vegan::vegdist(m, method = "bray"))[1, 2]
    data.frame(name1 = paste0(y, "_warming"),
               name2 = paste0(y, "_unwarming"), bray = d,
               stringsAsFactors = FALSE)
  }))
}

param_key <- c("Immigration", "Interaction", "Growth", "Noise")
md <- read_sample_metadata()

for (group in GROUPS) {
  if (!group %in% cfg$Group) stop("Unknown group: ", group)
  counts_emp <- read_empirical_counts(group, md)

  empirical_temporal <- list(
    warming = empirical_temporal_target(counts_emp, md, "warming"),
    unwarming = empirical_temporal_target(counts_emp, md, "unwarming")
  )
  empirical_spatial <- empirical_spatial_target(counts_emp, md)
  empirical_ratio <- empirical_richness_ratio_target(counts_emp, md)

  target_dir <- file.path(ROOT, "Calibration", "empirical_targets", group)
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(empirical_temporal$warming,
            file.path(target_dir, "temporal_bray_warming.csv"), row.names = FALSE)
  write.csv(empirical_temporal$unwarming,
            file.path(target_dir, "temporal_bray_unwarming.csv"), row.names = FALSE)
  write.csv(empirical_spatial,
            file.path(target_dir, "spatial_bray_warming_vs_unwarming.csv"), row.names = FALSE)
  write.csv(data.frame(Year = names(empirical_ratio), Ratio = as.numeric(empirical_ratio)),
            file.path(target_dir, "richness_ratio_warming_vs_unwarming.csv"), row.names = FALSE)

  meta_by_condition <- list()
  for (condition in c("warming", "unwarming")) {
    idir <- file.path(ROOT, "Resample", group, "raw", condition)
    ff <- sort(list.files(idir,
      pattern = "^Time_series_from_SLM_immigration_.*\\.csv$", full.names = TRUE))
    if (!length(ff)) stop("No calibration resamples in: ", idir)
    mm <- do.call(rbind, lapply(ff, parse_calibration_filename))
    mm$File <- ff
    meta_by_condition[[condition]] <- mm
  }

  match_cols <- c(param_key, "Seed")
  paired <- merge(meta_by_condition$warming, meta_by_condition$unwarming,
                  by = match_cols, suffixes = c("_warming", "_unwarming"))
  if (!nrow(paired)) stop(group, ": no matched warmed/unwarmed calibration files.")

  rows <- vector("list", nrow(paired))
  for (i in seq_len(nrow(paired))) {
    w <- read_counts(paired$File_warming[i])
    u <- read_counts(paired$File_unwarming[i])

    temporal_w <- mse_pair_tables(empirical_temporal$warming, bray_pairs(w),
                                  paste(group, "warming temporal"))
    temporal_u <- mse_pair_tables(empirical_temporal$unwarming, bray_pairs(u),
                                  paste(group, "unwarming temporal"))
    spatial <- mse_pair_tables(empirical_spatial, spatial_pairs(w, u),
                               paste(group, "warming-vs-unwarming spatial"))

    rw <- colSums(w > 0)
    ru <- colSums(u > 0)
    ratio <- rw / ru
    nratio <- min(length(ratio), length(empirical_ratio))
    richness_mse <- mean((empirical_ratio[seq_len(nratio)] - ratio[seq_len(nratio)])^2,
                         na.rm = TRUE)

    rows[[i]] <- data.frame(
      Group = group,
      paired[i, match_cols, drop = FALSE],
      TemporalBray_warming = temporal_w,
      TemporalBray_unwarming = temporal_u,
      TemporalBrayDifference = temporal_w + temporal_u,
      SpatialBrayDifference = spatial,
      RichnessRatioDifference = richness_mse,
      stringsAsFactors = FALSE
    )
  }

  seed_tab <- do.call(rbind, rows)
  outdir <- file.path(ROOT, "Calibration", "community")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  write.csv(seed_tab, file.path(outdir, paste0(group, "_community_metrics_seed.csv")),
            row.names = FALSE, quote = FALSE)

  aggregate_formula <- as.formula(
    paste("cbind(TemporalBray_warming, TemporalBray_unwarming,",
          "TemporalBrayDifference, SpatialBrayDifference, RichnessRatioDifference) ~",
          paste(param_key, collapse = " + "))
  )
  param_tab <- aggregate(aggregate_formula, data = seed_tab, FUN = mean, na.rm = TRUE)
  param_tab$Group <- group
  param_tab <- param_tab[, c("Group", param_key,
                             "TemporalBray_warming", "TemporalBray_unwarming",
                             "TemporalBrayDifference", "SpatialBrayDifference",
                             "RichnessRatioDifference")]
  write.csv(param_tab, file.path(outdir, paste0(group, "_community_metrics_parameter.csv")),
            row.names = FALSE, quote = FALSE)

  message("Finished community calibration metrics: ", group)
}

message("All sGLV community calibration metrics finished.")
