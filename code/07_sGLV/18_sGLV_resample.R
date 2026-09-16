# ==============================================================================
# 18_sGLV_resample.R
# Resample simulated relative-abundance trajectories to group-specific library size
# for both calibration-grid and top-10 process-intensity simulations.
#
# Environment variables:
#   TPL_PROJECT_ROOT     repository root (default: getwd())
#   TPL_SGLV_DATA_DIR    model-derived input root (default: ./data/sGLV)
#   TPL_SGLV_RESAMPLE_MODE  calibration | process_intensity | all (default: all)
#   TPL_SGLV_TPL         tTTPL | pTTPL | both (process-intensity mode; default: both)
#   TPL_SGLV_GROUPS      comma-separated groups or all (default: all)
#   TPL_SGLV_RESAMPLE_SEED integer base seed (default: 20260915)
#   TPL_SGLV_OVERWRITE   TRUE/FALSE (default: FALSE)
#
# Calibration:
#   get_timeseries/<Group>/raw/<condition>/ -> Resample/<Group>/raw/<condition>/
#
# Process intensity:
#   get_timeseries/<Group>/<TPL>/process_intensity_top10/<process>/<condition>/
#   -> Resample/<Group>/<TPL>/process_intensity_top10/<process>/<condition>/
# ==============================================================================

rm(list = ls())

PROJECT_DIR <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(Sys.getenv("TPL_SGLV_WORK_DIR",
                                 unset = file.path(PROJECT_DIR, "results", "17_sGLV")),
                      winslash = "/", mustWork = FALSE)
DATA_SGLV <- normalizePath(Sys.getenv("TPL_SGLV_DATA_DIR",
                                      unset = file.path(PROJECT_DIR, "data", "sGLV")),
                           winslash = "/", mustWork = FALSE)
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

SCRIPT_DIR <- file.path(PROJECT_DIR, "code", "07_sGLV")
CONFIG_FILE <- file.path(SCRIPT_DIR, "sGLV_group_config.csv")
if (!file.exists(CONFIG_FILE)) stop("Cannot find sGLV_group_config.csv")

cfg <- read.csv(CONFIG_FILE, check.names = FALSE, stringsAsFactors = FALSE)

parse_groups <- function(x) {
  if (!nzchar(x) || tolower(x) == "all") return(cfg$Group)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

GROUPS <- parse_groups(Sys.getenv("TPL_SGLV_GROUPS", unset = "all"))
MODE <- Sys.getenv("TPL_SGLV_RESAMPLE_MODE", unset = "all")
TPL_MODE <- Sys.getenv("TPL_SGLV_TPL", unset = "both")
BASE_SEED <- as.integer(Sys.getenv("TPL_SGLV_RESAMPLE_SEED", unset = "20260915"))
OVERWRITE <- toupper(Sys.getenv("TPL_SGLV_OVERWRITE", unset = "FALSE")) == "TRUE"

if (!MODE %in% c("calibration", "process_intensity", "all"))
  stop("TPL_SGLV_RESAMPLE_MODE must be calibration, process_intensity, or all")
if (!TPL_MODE %in% c("tTTPL", "pTTPL", "both"))
  stop("TPL_SGLV_TPL must be tTTPL, pTTPL, or both")

conditions <- c("warming", "unwarming")
processes <- c("immigration", "growth", "interaction", "noise")
tpls <- if (TPL_MODE == "both") c("tTTPL", "pTTPL") else TPL_MODE

read_taxa <- function(group, condition, depth) {
  f <- file.path(
    DATA_SGLV, "rawdata",
    paste0("resample", depth, "_otutab_", group, "_", condition, "_2010_mean.csv")
  )
  if (!file.exists(f)) stop("Missing empirical initial-abundance file: ", f)
  x <- read.csv(f, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(x)
}

read_simulation <- function(file, expected_taxa) {
  x <- read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  # Canonical Python output: years x taxa. Legacy outputs may be taxa x time.
  if (all(expected_taxa %in% colnames(x))) {
    x <- x[, expected_taxa, drop = FALSE]
    mat <- t(as.matrix(x))
  } else if (all(expected_taxa %in% rownames(x))) {
    x <- x[expected_taxa, , drop = FALSE]
    mat <- as.matrix(x)
  } else if (ncol(x) == length(expected_taxa)) {
    colnames(x) <- expected_taxa
    mat <- t(as.matrix(x))
  } else if (nrow(x) == length(expected_taxa)) {
    rownames(x) <- expected_taxa
    mat <- as.matrix(x)
  } else {
    stop("Cannot align simulated taxa in: ", file)
  }
  storage.mode(mat) <- "double"
  mat[!is.finite(mat) | mat < 0] <- 0
  rownames(mat) <- expected_taxa
  if (ncol(mat) < 11) stop("Fewer than 11 simulated snapshots in: ", file)
  colnames(mat) <- paste0("Y", 2010 + seq_len(ncol(mat)) - 1L)
  mat
}

resample_one <- function(infile, outfile, taxa, depth, seed) {
  if (file.exists(outfile) && !OVERWRITE) return(FALSE)
  mat <- read_simulation(infile, taxa)
  set.seed(seed)
  out <- matrix(0L, nrow = length(taxa), ncol = ncol(mat),
                dimnames = list(taxa, colnames(mat)))
  for (j in seq_len(ncol(mat))) {
    p <- mat[, j]
    if (!is.finite(sum(p)) || sum(p) <= 0) {
      stop("Non-positive simulated total abundance in ", infile, " at ", colnames(mat)[j])
    }
    draw <- sample.int(length(taxa), size = depth, replace = TRUE, prob = p)
    out[, j] <- tabulate(draw, nbins = length(taxa))
  }
  dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(OTUs = rownames(out), out, check.names = FALSE),
            outfile, row.names = FALSE, quote = FALSE)
  TRUE
}

run_directory <- function(indir, outdir, group, condition, depth, pattern) {
  if (!dir.exists(indir)) {
    message("Skip missing directory: ", indir)
    return(invisible(0L))
  }
  files <- sort(list.files(indir, pattern = pattern, full.names = TRUE))
  if (!length(files)) {
    message("No matching simulations in: ", indir)
    return(invisible(0L))
  }
  taxa <- read_taxa(group, condition, depth)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  nw <- 0L
  for (i in seq_along(files)) {
    outfile <- file.path(outdir, basename(files[i]))
    # Stable file-specific seed so reruns are reproducible and independent of run order.
    this_seed <- BASE_SEED + match(group, cfg$Group) * 100000L +
      match(condition, conditions) * 10000L + i
    nw <- nw + as.integer(resample_one(files[i], outfile, taxa, depth, this_seed))
  }
  message(group, " | ", condition, " | ", basename(indir), ": wrote ", nw,
          " of ", length(files), " files")
  invisible(nw)
}

for (group in GROUPS) {
  ii <- match(group, cfg$Group)
  if (is.na(ii)) stop("Unknown group: ", group)
  depth <- as.integer(cfg$ResampleDepth[ii])

  if (MODE %in% c("calibration", "all")) {
    for (condition in conditions) {
      run_directory(
        file.path(ROOT, "get_timeseries", group, "raw", condition),
        file.path(ROOT, "Resample", group, "raw", condition),
        group, condition, depth,
        "^Time_series_from_SLM_immigration_.*\\.csv$"
      )
    }
  }

  if (MODE %in% c("process_intensity", "all")) {
    for (tpl in tpls) {
      for (process in processes) {
        for (condition in conditions) {
          run_directory(
            file.path(ROOT, "get_timeseries", group, tpl,
                      "process_intensity_top10", process, condition),
            file.path(ROOT, "Resample", group, tpl,
                      "process_intensity_top10", process, condition),
            group, condition, depth,
            "^Time_series_from_SLM_top10rank.*\\.csv$"
          )
        }
      }
    }
  }
}

message("sGLV resampling finished.")
