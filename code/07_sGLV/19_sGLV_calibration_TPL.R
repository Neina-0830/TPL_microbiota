# ==============================================================================
# 19_sGLV_calibration_TPL.R
# Recalculate tTTPL and pTTPL for all calibration-grid simulations.
#
# Requires resampled files from 18_sGLV_resample.R.
#
# Run from the repository root. External inputs are under ./data;
# generated phylogenetic distances are read from ./results/pairwise_distance.
# Optional environment variables:
#   TPL_PROJECT_ROOT
#   TPL_SGLV_GROUPS = all or comma-separated
#   TPL_SGLV_CALIBRATION_TPL = tTTPL | pTTPL | both (default: both)
# ==============================================================================

rm(list = ls())

PROJECT_ROOT <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()),
                              winslash = "/", mustWork = TRUE)
ROOT <- normalizePath(Sys.getenv("TPL_SGLV_WORK_DIR",
                                 unset = file.path(PROJECT_ROOT, "results", "17_sGLV")),
                      winslash = "/", mustWork = FALSE)
PHYLO_ROOT <- file.path(PROJECT_ROOT, "results", "pairwise_distance")
dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(ROOT)

needed <- c("iCAMP", "bigmemory")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install required R packages: ", paste(miss, collapse = ", "))

SCRIPT_DIR <- file.path(PROJECT_ROOT, "code", "07_sGLV")
CONFIG_FILE <- file.path(SCRIPT_DIR, "sGLV_group_config.csv")
cfg <- read.csv(CONFIG_FILE, check.names = FALSE, stringsAsFactors = FALSE)

GROUPS <- Sys.getenv("TPL_SGLV_GROUPS", unset = "all")
GROUPS <- if (tolower(GROUPS) == "all") cfg$Group else trimws(strsplit(GROUPS, ",")[[1]])
TPL_MODE <- Sys.getenv("TPL_SGLV_CALIBRATION_TPL", unset = "both")
if (!TPL_MODE %in% c("tTTPL", "pTTPL", "both")) stop("Invalid TPL_SGLV_CALIBRATION_TPL")
TPLS <- if (TPL_MODE == "both") c("tTTPL", "pTTPL") else TPL_MODE

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

read_resampled <- function(file) {
  x <- read.csv(file, row.names = 1, check.names = FALSE)
  # 18_sGLV_resample.R writes OTUs as the first column -> row.names=1 gives taxa x years.
  if (ncol(x) < 11) stop("Fewer than 11 annual columns in: ", file)
  x <- x[, seq_len(11), drop = FALSE]
  x <- as.matrix(t(x))  # years x taxa
  storage.mode(x) <- "double"
  x
}

fit_tpl <- function(mat) {
  mu <- colMeans(mat)
  va <- apply(mat, 2, var)
  keep <- is.finite(mu) & is.finite(va) & mu > 0 & va > 0
  if (sum(keep) < 3) return(data.frame(
    Intercept = NA_real_, Beta = NA_real_, R2 = NA_real_, P = NA_real_,
    N_taxa = sum(keep)
  ))
  dd <- data.frame(LogMean = log10(mu[keep]), LogVar = log10(va[keep]))
  fit <- lm(LogVar ~ LogMean, data = dd)
  ss <- summary(fit)
  data.frame(
    Intercept = unname(coef(fit)[1]),
    Beta = unname(coef(fit)[2]),
    R2 = unname(ss$r.squared),
    P = unname(ss$coefficients[2, 4]),
    N_taxa = nrow(dd)
  )
}

mpdi.big <- function(comm, pd.desc = "pd.desc", pd.spname = NULL,
                     pd.wd = getwd(), sp.limit = 10000,
                     spname.check = TRUE, abundance.weighted = TRUE) {
  if (spname.check) {
    if (is.null(pd.spname)) {
      pd.spname <- read.csv(file.path(pd.wd, "pd.taxon.name.csv"),
                            row.names = 1, stringsAsFactors = FALSE)[, 1]
    }
    checked <- iCAMP::match.name(name.check = pd.spname, cn.list = list(comm = comm))
    comm <- checked$comm
  }
  ids <- match(colnames(comm), pd.spname)
  if (anyNA(ids)) stop("Some taxa could not be matched to the phylogenetic-distance matrix.")
  if (!abundance.weighted) comm[comm > 0] <- 1
  rs <- rowSums(comm)
  if (any(!is.finite(rs)) || any(rs <= 0)) stop("Non-positive community total.")
  comm <- as.matrix(comm / rs)

  desc <- file.path(pd.wd, pd.desc)
  pd <- try(bigmemory::attach.big.matrix(dget(desc)), silent = TRUE)
  if (inherits(pd, "try-error")) pd <- bigmemory::attach.big.matrix(desc)

  if (nrow(pd) <= sp.limit) {
    comd <- comm %*% pd[ids, ids]
  } else {
    nsp <- ncol(comm)
    chunk <- max(1L, floor(sp.limit * sp.limit / nsp))
    starts <- seq(1L, nsp, by = chunk)
    comd <- matrix(0, nrow(comm), nsp)
    for (from in starts) {
      to <- min(from + chunk - 1L, nsp)
      comd[, from:to] <- comm %*% pd[ids, ids[from:to], drop = FALSE]
    }
  }
  ans <- comd * comm
  colnames(ans) <- colnames(comm)
  ans
}

pdist_candidates <- function(group) {
  phy_group <- if (group == "Fungi_18S") "Fungi" else group
  file.path(PHYLO_ROOT, phy_group, "All_OTUs_2010to2020_pdist_out.rda")
}

load_pdist <- function(group) {
  cand <- pdist_candidates(group)
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop(
    "No pairwise-distance RDA found for ", group, ". Tried:\n",
    paste(cand, collapse = "\n")
  )
  e <- new.env(parent = emptyenv())
  load(hit[1], envir = e)
  if (!exists("pdist.out", envir = e, inherits = FALSE))
    stop("Object pdist.out not found in: ", hit[1])
  obj <- get("pdist.out", envir = e)
  if (!file.exists(file.path(obj$pd.wd, obj$pd.file))) {
    local_desc <- file.path(dirname(hit[1]), obj$pd.file)
    if (!file.exists(local_desc)) stop("Cannot locate pd.desc for ", group)
    obj$pd.wd <- dirname(hit[1])
  }
  obj
}

calc_pTTPL <- function(mat, pdist.out) {
  community_mpd <- iCAMP::mpd.big(
    mat, pd.desc = pdist.out$pd.file, pd.spname = pdist.out$tip.label,
    pd.wd = pdist.out$pd.wd, spname.check = TRUE
  )
  contrib <- mpdi.big(
    mat, pd.desc = pdist.out$pd.file, pd.spname = pdist.out$tip.label,
    pd.wd = pdist.out$pd.wd, spname.check = TRUE
  )
  mpdik <- contrib / community_mpd
  if (any(!is.finite(mpdik))) stop("Non-finite MPD_i,k generated.")
  list(beta = fit_tpl(mpdik), mpdik = mpdik)
}

for (group in GROUPS) {
  if (!group %in% cfg$Group) stop("Unknown group: ", group)
  pdist.out <- if ("pTTPL" %in% TPLS) load_pdist(group) else NULL

  for (condition in c("warming", "unwarming")) {
    indir <- file.path(ROOT, "Resample", group, "raw", condition)
    files <- sort(list.files(indir,
      pattern = "^Time_series_from_SLM_immigration_.*\\.csv$", full.names = TRUE))
    if (!length(files)) stop("No calibration resample files found in: ", indir)

    t_rows <- list()
    p_rows <- list()
    for (k in seq_along(files)) {
      meta <- parse_calibration_filename(files[k])
      mat <- read_resampled(files[k])

      if ("tTTPL" %in% TPLS) {
        rr <- cbind(meta, fit_tpl(mat), Condition = condition, Group = group)
        t_rows[[length(t_rows) + 1L]] <- rr
      }

      if ("pTTPL" %in% TPLS) {
        pp <- calc_pTTPL(mat, pdist.out)
        rr <- cbind(meta, pp$beta, Condition = condition, Group = group)
        p_rows[[length(p_rows) + 1L]] <- rr
      }
    }

    outdir <- file.path(ROOT, "Taylor", group, "raw", condition)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    if (length(t_rows)) {
      write.csv(do.call(rbind, t_rows),
        file.path(outdir,
          "SLM_all_immigration_interaction_growth_all_seed_noise_all_species_Taylor_beta.csv"),
        row.names = FALSE, quote = FALSE)
    }
    if (length(p_rows)) {
      write.csv(do.call(rbind, p_rows),
        file.path(outdir,
          "MPD_weighted_SLM_all_immigration_interaction_growth_all_seed_noise_all_species_Taylor_beta.csv"),
        row.names = FALSE, quote = FALSE)
    }
    message("Finished calibration TPL: ", group, " | ", condition)
  }
}

message("Calibration tTTPL/pTTPL fitting finished.")
