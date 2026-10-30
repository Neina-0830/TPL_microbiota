# ==============================================================================
# iCAMP community-assembly analysis and process-summary workflow
# ==============================================================================
# This cleaned script is based on the supplied Bacteria iCAMP workflow and now
# includes the raw-data -> iCAMP -> process-summary stages that were absent from
# the earlier cleaned release.
#
# The supplied source explicitly used:
#   rand.time       = 200
#   bin.size.limit  = 24
#   ds              = 0.2
#   sig.index       = "Confidence"
#   phylo.metric    = "bMPD"
#   phylo.rand.scale = "within.bin"
#   taxa.rand.scale  = "across.all"
#   ses.cut         = 1.96
#   rc.cut          = 0.95
#   conf.cut        = 0.975
#   omit.option     = "no"
#
# The optional niche-difference / phylogenetic-signal exploration (historical
# steps 7-8) remains omitted because it was commented out in the supplied source;
# the selected bin size used in the actual run is preserved directly above.
#
# By default this script is configured for the supplied Bacteria analysis. The
# same code can be reused for another soil-biotic group by overriding the input
# files and GROUP/PREFIX through environment variables, but those group-specific
# filenames/settings must come from the corresponding original analysis and are
# not inferred here.
#
# Two execution modes are supported:
#   1) Full run from raw community/tree data (default when no detail object exists)
#   2) Reuse an existing iCAMP Confidence detail object and rerun only summaries
#
# Environment variables (all optional):
#   TPL_PROJECT_ROOT        repository/project root; default = getwd()
#   TPL_ICAMP_GROUP         default = Bacteria
#   TPL_ICAMP_DIR           default = <root>/iCAMP/<Group>
#   TPL_ICAMP_OUTPUT_DIR    default = <iCAMP dir>/Outputs
#   TPL_ICAMP_COMM          OTU/ASV table
#   TPL_ICAMP_TREE          rooted/unrooted Newick tree
#   TPL_ICAMP_CLASS         classification table
#   TPL_ICAMP_TREAT         treatment table
#   TPL_ICAMP_ENV           environmental table; optional for the active workflow
#   TPL_ICAMP_DETAIL        saved *.iCAMP.Confidence.detail.rda to reuse
#   TPL_ICAMP_RUN_RAW       TRUE/FALSE; default = TRUE only when no detail exists
#   TPL_ICAMP_RAND_TIME     default = 200
#   TPL_ICAMP_NWORKER       default = 20
#   TPL_ICAMP_MEMORY_G      default = 50
#   TPL_ICAMP_TREAT_COLUMN  1-based column used by icamp.boot; default = 1
#
# Main outputs:
#   <Group>.iCAMP.Confidence.detail.cleaned.rds
#   <Group>.iCAMP.Summary_Confidence.rda
#   <Group>.ProcessImportance_EachGroup_Confidence.csv
#   <Group>.ProcessImportance_EachBin_EachGroup_Confidence.csv
#   <Group>.ProcessImportance_EachTurnover_Confidence.csv
#   <Group>.BinContributeToProcess_EachGroup_Confidence.csv
#   <Group>.Taxon_Bin_Confidence.csv
#   <Group>.Bin_TopTaxon_Confidence.csv
#   <Group>.iCAMP.Boot*.rda/csv
#   <Group>.iCAMP_settings_used.csv
# ============================================================================== 

rm(list = ls())
t0 <- Sys.time()

# ------------------------------------------------------------------------------
# 0. Helpers
# ------------------------------------------------------------------------------

env_bool <- function(name, default = FALSE) {
  x <- Sys.getenv(name, unset = "")
  if (!nzchar(x)) return(default)
  x <- tolower(trimws(x))
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop(name, " must be TRUE/FALSE (or 1/0).")
}

as_int_env <- function(name, default) {
  x <- Sys.getenv(name, unset = as.character(default))
  out <- suppressWarnings(as.integer(x))
  if (!is.finite(out)) stop(name, " must be an integer.")
  out
}

as_num_env <- function(name, default) {
  x <- Sys.getenv(name, unset = as.character(default))
  out <- suppressWarnings(as.numeric(x))
  if (!is.finite(out)) stop(name, " must be numeric.")
  out
}

assert_file <- function(file, label) {
  if (!file.exists(file)) stop(label, " not found: ", file)
}

read_table <- function(file) {
  utils::read.csv(
    file,
    header = TRUE,
    sep = ",",
    row.names = 1,
    as.is = TRUE,
    stringsAsFactors = FALSE,
    comment.char = "",
    check.names = FALSE
  )
}

load_icamp_detail_rda <- function(file) {
  e <- new.env(parent = emptyenv())
  loaded <- load(file, envir = e)
  objs <- lapply(loaded, function(nm) get(nm, envir = e, inherits = FALSE))
  ok <- vapply(
    objs,
    function(x) is.list(x) && !is.null(x$detail),
    logical(1)
  )
  if (!any(ok)) {
    stop("No iCAMP object with a $detail component was found in: ", file)
  }
  objs[[which(ok)[1]]]
}

# ------------------------------------------------------------------------------
# 1. Project and input configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv("TPL_PROJECT_ROOT", unset = getwd())
if (!dir.exists(PROJECT_DIR)) stop("TPL_PROJECT_ROOT does not exist: ", PROJECT_DIR)
PROJECT_DIR <- normalizePath(PROJECT_DIR, mustWork = TRUE)
setwd(PROJECT_DIR)

GROUP <- Sys.getenv("TPL_ICAMP_GROUP", unset = "Bacteria")
PREFIX <- GROUP

DATA_DIR <- file.path(PROJECT_DIR, "data")
PREP_DIR <- file.path(DATA_DIR, "prepared")
OUT_DIR <- Sys.getenv(
  "TPL_ICAMP_OUTPUT_DIR",
  unset = file.path(PROJECT_DIR, "results", "14_iCAMP", GROUP)
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COMM_FILES <- c(
  Archaea  = "resample7842_otutab_Archaea_2010to2020.csv",
  Bacteria = "resample26890_otutab_Bacteria_2010to2020.csv",
  Fungi    = "resample541_otutab_Fungi_2010to2020.csv",
  Protists = "resample1064_otutab_Protists_2010to2020.csv",
  Nematoda = "resample296_otutab_Nematoda_2010to2020.csv"
)
TREE_FILES <- c(
  Archaea  = "rooted_tree_otus_Archaea_2010to2020_resample7842.nwk",
  Bacteria = "rooted_tree_otus_Bacteria_2010to2020_resample26890.nwk",
  Fungi    = "rooted_tree_otus_Fungi_2010to2020_resample541.nwk",
  Protists = "rooted_tree_otus_Protists_2010to2020_resample1064.nwk",
  Nematoda = "rooted_tree_otus_Nematoda_2010to2020_resample296.nwk"
)
CLASS_FILES <- c(
  Archaea  = "Archaea_taxonomy.csv",
  Bacteria = "Bacteria_taxonomy.csv",
  Fungi    = "Fungi_taxonomy.csv",
  Protists = "Protists_taxonomy.csv",
  Nematoda = "Nematoda_taxonomy.csv"
)

if (!GROUP %in% names(COMM_FILES)) {
  stop("TPL_ICAMP_GROUP must be one of: ", paste(names(COMM_FILES), collapse = ", "))
}

COMM_FILE <- Sys.getenv(
  "TPL_ICAMP_COMM",
  unset = file.path(DATA_DIR, unname(COMM_FILES[GROUP]))
)
TREE_FILE <- Sys.getenv(
  "TPL_ICAMP_TREE",
  unset = file.path(DATA_DIR, unname(TREE_FILES[GROUP]))
)
CLASS_FILE <- Sys.getenv(
  "TPL_ICAMP_CLASS",
  unset = file.path(PREP_DIR, "taxonomy", unname(CLASS_FILES[GROUP]))
)
TREAT_FILE <- Sys.getenv(
  "TPL_ICAMP_TREAT",
  unset = file.path(PREP_DIR, "iCAMP_treatment_by_sample.csv")
)
ENV_FILE <- Sys.getenv(
  "TPL_ICAMP_ENV",
  unset = file.path(PREP_DIR, "environment_metadata.csv")
)

DETAIL_FILE <- Sys.getenv(
  "TPL_ICAMP_DETAIL",
  unset = file.path(OUT_DIR, paste0(PREFIX, ".iCAMP.Confidence.detail.rda"))
)
CLEAN_DETAIL_RDS <- file.path(
  OUT_DIR,
  paste0(PREFIX, ".iCAMP.Confidence.detail.cleaned.rds")
)

RAND_TIME <- as_int_env("TPL_ICAMP_RAND_TIME", 200L)
NWORKER <- as_int_env("TPL_ICAMP_NWORKER", 20L)
MEMORY_G <- as_num_env("TPL_ICAMP_MEMORY_G", 50)
# The prepared treatment table has Treatment as its first data column.
TREAT_COLUMN <- as_int_env("TPL_ICAMP_TREAT_COLUMN", 1L)

# Scientific settings supplied in the original analysis.
BIN_SIZE_LIMIT <- 24L
DS <- 0.2
SIG_INDEX <- "Confidence"
PHYLO_METRIC <- "bMPD"
PHYLO_RAND_SCALE <- "within.bin"
TAXA_RAND_SCALE <- "across.all"
SES_CUT <- 1.96
RC_CUT <- 0.95
CONF_CUT <- 0.975
OMIT_OPTION <- "no"

DETAIL_EXISTS <- file.exists(CLEAN_DETAIL_RDS) || file.exists(DETAIL_FILE)
RUN_RAW <- env_bool("TPL_ICAMP_RUN_RAW", default = !DETAIL_EXISTS)

# ------------------------------------------------------------------------------
# 2. Packages
# ------------------------------------------------------------------------------

required_pkgs <- c("iCAMP", "ape")
missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_pkgs)) {
  stop(
    "Install required R package(s) before running: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 3. Read treatment/classification tables (needed in both execution modes)
# ------------------------------------------------------------------------------

assert_file(TREAT_FILE, "Treatment table")
assert_file(CLASS_FILE, "Classification table")

treat <- read_table(TREAT_FILE)
clas <- read_table(CLASS_FILE)

if (TREAT_COLUMN < 1L || TREAT_COLUMN > ncol(treat)) {
  stop("TPL_ICAMP_TREAT_COLUMN is outside the treatment-table column range.")
}

# ------------------------------------------------------------------------------
# 4. Full iCAMP run from raw data, or load an existing Confidence detail object
# ------------------------------------------------------------------------------

if (RUN_RAW) {

  assert_file(COMM_FILE, "Community table")
  if (!nzchar(TREE_FILE)) {
    stop(
      "No default tree filename is known for GROUP='", GROUP,
      "'. Set TPL_ICAMP_TREE explicitly."
    )
  }
  assert_file(TREE_FILE, "Phylogenetic tree")

  comm <- t(read_table(COMM_FILE))
  tree <- ape::read.tree(file = TREE_FILE)

  # Preserve the original sample-matching behavior when the environmental table
  # is available. The active selected-bin iCAMP run itself does not require env.
  if (file.exists(ENV_FILE)) {
    env <- read_table(ENV_FILE)
    sampid_check <- iCAMP::match.name(
      rn.list = list(comm = comm, treat = treat, env = env)
    )
    comm <- sampid_check$comm
    treat <- sampid_check$treat
    env <- sampid_check$env
  } else {
    message("Environmental table not found; matching community and treatment tables only: ", ENV_FILE)
    sampid_check <- iCAMP::match.name(
      rn.list = list(comm = comm, treat = treat)
    )
    comm <- sampid_check$comm
    treat <- sampid_check$treat
  }

  # Remove ghost taxa after sample matching.
  comm <- comm[, colSums(comm) > 0, drop = FALSE]

  # Match OTU/ASV IDs across community table, classification and tree.
  spid_check <- iCAMP::match.name(
    cn.list = list(comm = comm),
    rn.list = list(clas = clas),
    tree.list = list(tree = tree)
  )
  comm <- spid_check$comm
  clas <- spid_check$clas
  tree <- spid_check$tree

  # Pairwise phylogenetic distances. Reuse the big-memory files if already made.
  pd_desc <- file.path(OUT_DIR, "pd.desc")
  pd_taxon <- file.path(OUT_DIR, "pd.taxon.name.csv")

  if (!file.exists(pd_desc) || !file.exists(pd_taxon)) {
    pd.big <- iCAMP::pdist.big(
      tree = tree,
      wd = OUT_DIR,
      nworker = NWORKER,
      memory.G = MEMORY_G
    )
  } else {
    pd.big <- list(
      tip.label = utils::read.csv(
        pd_taxon,
        row.names = 1,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )[, 1],
      pd.wd = OUT_DIR,
      pd.file = "pd.desc",
      pd.name.file = "pd.taxon.name.csv"
    )
  }

  # The supplied source midpoint-rooted an unrooted tree before icamp.big.
  if (!ape::is.rooted(tree)) {
    tree.rt <- iCAMP::midpoint.root.big(
      tree = tree,
      pd.desc = pd.big$pd.file,
      pd.spname = pd.big$tip.label,
      pd.wd = pd.big$pd.wd,
      nworker = NWORKER
    )
    tree <- tree.rt$tree
  }

  message("Running iCAMP for ", GROUP, " ...")
  message(
    "Settings: rand.time=", RAND_TIME,
    ", bin.size.limit=", BIN_SIZE_LIMIT,
    ", ds=", DS,
    ", sig.index=", SIG_INDEX
  )

  icres <- iCAMP::icamp.big(
    comm = comm,
    pd.desc = pd.big$pd.file,
    pd.spname = pd.big$tip.label,
    pd.wd = pd.big$pd.wd,
    rand = RAND_TIME,
    tree = tree,
    prefix = PREFIX,
    ds = DS,
    pd.cut = NA,
    sp.check = TRUE,
    phylo.rand.scale = PHYLO_RAND_SCALE,
    taxa.rand.scale = TAXA_RAND_SCALE,
    phylo.metric = PHYLO_METRIC,
    sig.index = SIG_INDEX,
    bin.size.limit = BIN_SIZE_LIMIT,
    nworker = NWORKER,
    memory.G = MEMORY_G,
    rtree.save = FALSE,
    detail.save = TRUE,
    qp.save = FALSE,
    detail.null = FALSE,
    ignore.zero = TRUE,
    output.wd = OUT_DIR,
    correct.special = TRUE,
    unit.sum = rowSums(comm),
    special.method = "depend",
    ses.cut = SES_CUT,
    rc.cut = RC_CUT,
    conf.cut = CONF_CUT,
    omit.option = OMIT_OPTION,
    meta.ab = NULL
  )

  # Store a canonical object with a stable RDS name. iCAMP itself may also write
  # its standard *.iCAMP.Confidence.detail.rda file into OUT_DIR.
  saveRDS(icres, CLEAN_DETAIL_RDS)

} else {

  if (file.exists(CLEAN_DETAIL_RDS)) {
    message("Reusing cleaned iCAMP detail object: ", CLEAN_DETAIL_RDS)
    icres <- readRDS(CLEAN_DETAIL_RDS)
  } else {
    assert_file(DETAIL_FILE, "iCAMP Confidence detail RDA")
    message("Reusing saved iCAMP detail object: ", DETAIL_FILE)
    icres <- load_icamp_detail_rda(DETAIL_FILE)
    saveRDS(icres, CLEAN_DETAIL_RDS)
  }
}

if (!is.list(icres) || is.null(icres$detail)) {
  stop("The iCAMP result object does not contain the required $detail component.")
}

# ------------------------------------------------------------------------------
# 5. iCAMP bin/process summaries
# ------------------------------------------------------------------------------

icbin <- iCAMP::icamp.bins(
  icamp.detail = icres$detail,
  treat = treat,
  clas = clas,
  silent = FALSE,
  boot = TRUE,
  rand.time = RAND_TIME,
  between.group = TRUE
)

save(
  icbin,
  file = file.path(OUT_DIR, paste0(PREFIX, ".iCAMP.Summary_Confidence.rda"))
)

utils::write.csv(
  icbin$Pt,
  file.path(OUT_DIR, paste0(PREFIX, ".ProcessImportance_EachGroup_Confidence.csv")),
  row.names = FALSE
)
utils::write.csv(
  icbin$Ptk,
  file.path(OUT_DIR, paste0(PREFIX, ".ProcessImportance_EachBin_EachGroup_Confidence.csv")),
  row.names = FALSE
)
utils::write.csv(
  icbin$Ptuv,
  file.path(OUT_DIR, paste0(PREFIX, ".ProcessImportance_EachTurnover_Confidence.csv")),
  row.names = FALSE
)
utils::write.csv(
  icbin$BPtk,
  file.path(OUT_DIR, paste0(PREFIX, ".BinContributeToProcess_EachGroup_Confidence.csv")),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    ID = rownames(icbin$Class.Bin),
    icbin$Class.Bin,
    stringsAsFactors = FALSE
  ),
  file.path(OUT_DIR, paste0(PREFIX, ".Taxon_Bin_Confidence.csv")),
  row.names = FALSE
)
utils::write.csv(
  icbin$Bin.TopClass,
  file.path(OUT_DIR, paste0(PREFIX, ".Bin_TopTaxon_Confidence.csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Bootstrapping test for the selected treatment column
# ------------------------------------------------------------------------------

if (is.null(icres$CbMPDiCBraya)) {
  stop("The iCAMP result object does not contain $CbMPDiCBraya required by icamp.boot.")
}

treat.use <- treat[, TREAT_COLUMN, drop = FALSE]
tr_name <- colnames(treat)[TREAT_COLUMN]

icboot <- iCAMP::icamp.boot(
  icamp.result = icres$CbMPDiCBraya,
  treat = treat.use,
  rand.time = RAND_TIME,
  compare = TRUE,
  silent = FALSE,
  between.group = TRUE,
  ST.estimation = TRUE
)

save(
  icboot,
  file = file.path(OUT_DIR, paste0(PREFIX, ".iCAMP.Boot.", tr_name, ".rda"))
)
utils::write.csv(
  icboot$summary,
  file.path(OUT_DIR, paste0(PREFIX, ".iCAMP.BootSummary.", tr_name, ".csv")),
  row.names = FALSE
)
utils::write.csv(
  icboot$compare,
  file.path(OUT_DIR, paste0(PREFIX, ".iCAMP.Compare.", tr_name, ".csv")),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Save an explicit settings/provenance table
# ------------------------------------------------------------------------------

settings_used <- data.frame(
  Group = GROUP,
  Prefix = PREFIX,
  Run_raw_icamp = RUN_RAW,
  Rand_time = RAND_TIME,
  Bin_size_limit = BIN_SIZE_LIMIT,
  ds = DS,
  Sig_index = SIG_INDEX,
  Phylo_metric = PHYLO_METRIC,
  Phylo_rand_scale = PHYLO_RAND_SCALE,
  Taxa_rand_scale = TAXA_RAND_SCALE,
  SES_cut = SES_CUT,
  RC_cut = RC_CUT,
  Confidence_cut = CONF_CUT,
  Omit_option = OMIT_OPTION,
  Nworker = NWORKER,
  Memory_G = MEMORY_G,
  Treat_column = TREAT_COLUMN,
  Treat_column_name = tr_name,
  stringsAsFactors = FALSE
)

utils::write.csv(
  settings_used,
  file.path(OUT_DIR, paste0(PREFIX, ".iCAMP_settings_used.csv")),
  row.names = FALSE
)

elapsed <- difftime(Sys.time(), t0, units = "mins")
message("iCAMP workflow completed for ", GROUP, ".")
message("Output directory: ", OUT_DIR)
message("Elapsed time (min): ", round(as.numeric(elapsed), 2))
