# ==============================================================================
# Prepare plot-level MPD_i,k matrices for pTTPL
# ==============================================================================
# This script consolidates the original Calculate_MPDik.R workflow and writes
# the exact plot-level file structure consumed by 04_pTTPL_fit_full_precision.R.
#
# Inputs:
#   ./data/resample*_otutab_*_2010to2020.csv
#   ./data/rooted_tree_otus_*_2010to2020_resample*.nwk
#   ./data/prepared/environment_metadata.csv
#
# Outputs:
#   ./results/03_pTTPL_MPDik/<Group>/UPARSE/<Plot>_MPD_weighted_otu_All_abun.csv
#   ./results/pairwise_distance/<Group>/All_OTUs_2010to2020_pdist_out.rda
#
# The normalized taxon-specific contribution is
#   MPD_i,k = (taxon-specific abundance-weighted distance contribution) /
#             (community abundance-weighted MPD).
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
OUTPUT_ROOT <- file.path(PROJECT_DIR, "results", "03_pTTPL_MPDik")
PD_ROOT <- file.path(PROJECT_DIR, "results", "pairwise_distance")

required_packages <- c("iCAMP", "ape", "bigmemory")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R package(s): ", paste(missing_packages, collapse = ", "))
}

read_numeric_matrix <- function(file) {
  if (!file.exists(file)) stop("Input file not found: ", file)
  x <- read.csv(
    file, header = TRUE, row.names = 1, check.names = FALSE,
    stringsAsFactors = FALSE
  )
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) stop("Non-numeric/NA values in: ", file)
  x
}

mpdi_big <- function(comm, pdist.out, sp.limit = 10000L) {
  # comm: samples x taxa
  pd.spname <- pdist.out$tip.label
  ids <- match(colnames(comm), pd.spname)
  if (anyNA(ids)) {
    stop("Some community taxa could not be matched to the phylogenetic-distance matrix.")
  }

  rs <- rowSums(comm)
  if (any(!is.finite(rs)) || any(rs <= 0)) {
    stop("At least one sample has a non-positive/non-finite total abundance.")
  }
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

load_or_make_pdist <- function(tree, group) {
  pd_dir <- file.path(PD_ROOT, group)
  dir.create(pd_dir, recursive = TRUE, showWarnings = FALSE)
  rda_file <- file.path(pd_dir, "All_OTUs_2010to2020_pdist_out.rda")

  if (file.exists(rda_file)) {
    e <- new.env(parent = emptyenv())
    load(rda_file, envir = e)
    if (exists("pdist.out", envir = e, inherits = FALSE)) {
      obj <- get("pdist.out", envir = e, inherits = FALSE)
      desc_path <- file.path(obj$pd.wd, obj$pd.file)
      if (file.exists(desc_path)) return(obj)
    }
    message(group, ": existing pdist RDA is not portable/complete; regenerating.")
  }

  pdist.out <- iCAMP::pdist.big(tree = tree, wd = pd_dir)
  save(pdist.out, file = rda_file)
  pdist.out
}

find_input <- function(...) {
  candidates <- c(...)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("None of the candidate input files exists:\n", paste(candidates, collapse = "\n"))
  }
  hit[1]
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

if (!file.exists(META_FILE)) stop("Metadata file not found: ", META_FILE)
meta <- read.csv(META_FILE, check.names = FALSE, stringsAsFactors = FALSE)
needed_meta <- c("Sample_name", "plot", "year")
miss <- setdiff(needed_meta, names(meta))
if (length(miss)) stop("Metadata missing column(s): ", paste(miss, collapse = ", "))
meta$Sample_name <- as.character(meta$Sample_name)
meta$plot <- as.character(meta$plot)
meta$year <- as.numeric(meta$year)
meta <- meta[meta$year %in% 2010:2020, , drop = FALSE]

for (group in names(group_config)) {
  message("\n=== pTTPL MPD_i,k preparation: ", group, " ===")
  cfg <- group_config[[group]]

  comm_file <- find_input(file.path(PHYLO_DATA_DIR, cfg$comm))
  tree_file <- find_input(file.path(PHYLO_DATA_DIR, cfg$tree))

  comm <- read_numeric_matrix(comm_file)  # taxa x samples
  tree <- ape::read.tree(tree_file)

  common_taxa <- intersect(rownames(comm), tree$tip.label)
  if (length(common_taxa) < 3L) stop(group, ": fewer than 3 taxa shared by community and tree.")
  comm <- comm[rownames(comm) %in% common_taxa, , drop = FALSE]
  tree <- ape::keep.tip(tree, common_taxa)

  # Keep only samples represented in the 2010-2020 metadata.
  sample_ids <- intersect(colnames(comm), meta$Sample_name)
  if (!length(sample_ids)) stop(group, ": no community sample names match metadata Sample_name.")
  comm <- comm[, sample_ids, drop = FALSE]

  pdist.out <- load_or_make_pdist(tree, group)
  comm_samples <- t(comm)  # samples x taxa

  community_mpd <- iCAMP::mpd.big(
    comm_samples,
    pd.desc = pdist.out$pd.file,
    pd.spname = pdist.out$tip.label,
    pd.wd = pdist.out$pd.wd,
    spname.check = TRUE
  )
  each_taxon <- mpdi_big(comm_samples, pdist.out)
  mpdik <- sweep(each_taxon, 1L, community_mpd, "/")
  if (any(!is.finite(mpdik))) stop(group, ": non-finite MPD_i,k values produced.")

  out_dir <- file.path(OUTPUT_ROOT, group, "UPARSE")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plots <- unique(meta$plot)
  for (plot_id in plots) {
    md <- meta[meta$plot == plot_id, , drop = FALSE]
    md <- md[order(md$year), , drop = FALSE]
    ids <- md$Sample_name
    missing_samples <- setdiff(ids, rownames(mpdik))
    if (length(missing_samples)) {
      stop(group, " / ", plot_id, ": missing MPD_i,k sample(s): ",
           paste(missing_samples, collapse = ", "))
    }
    plot_mat <- t(mpdik[ids, , drop = FALSE])  # taxa x annual samples
    write.csv(
      plot_mat,
      file.path(out_dir, paste0(plot_id, "_MPD_weighted_otu_All_abun.csv")),
      quote = FALSE
    )
  }

  write.csv(
    mpdik,
    file.path(out_dir, paste0(group, "_MPDik_all_samples.csv")),
    quote = FALSE
  )

  message(group, ": wrote ", length(plots), " plot-level MPD_i,k matrices.")
}

message("\nAll pTTPL MPD_i,k plot-level inputs completed.")
