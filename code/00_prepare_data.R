# ==============================================================================
# Prepare repository-local inputs from the raw files in ./data
# ==============================================================================
# This script performs only file-format harmonization required by the downstream
# analyses. It does not alter the experimental design or statistical models.
#
# Raw inputs expected in ./data:
#   resample7842_otutab_Archaea_2010to2020.csv
#   resample26890_otutab_Bacteria_2010to2020.csv
#   resample541_otutab_Fungi_2010to2020.csv
#   resample1064_otutab_Protists_2010to2020.csv
#   resample296_otutab_Nematoda_2010to2020.csv
#   corresponding taxonomy files
#   rooted_tree_otus_*_2010to2020_resample*.nwk
#   summary of environmental data2010-2020.csv
#   treatment_list10_20.csv
#
# Generated inputs are written to ./data/prepared and are safe to delete and
# regenerate at any time.
# ==============================================================================

rm(list = ls())

PROJECT_DIR <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
DATA_DIR <- file.path(PROJECT_DIR, "data")
PREP_DIR <- file.path(DATA_DIR, "prepared")
PLOT_DIR <- file.path(PREP_DIR, "plot_abundance")
TAX_DIR <- file.path(PREP_DIR, "taxonomy")
FULL_DIR <- file.path(PREP_DIR, "full_relative_abundance")

dir.create(PREP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAX_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FULL_DIR, recursive = TRUE, showWarnings = FALSE)

ENV_RAW <- file.path(DATA_DIR, "summary of environmental data2010-2020.csv")
TREAT_RAW <- file.path(DATA_DIR, "treatment_list10_20.csv")

if (!file.exists(ENV_RAW)) stop("Missing input: ", ENV_RAW)
if (!file.exists(TREAT_RAW)) stop("Missing input: ", TREAT_RAW)

# ------------------------------------------------------------------------------
# 1. Harmonize sample-level metadata
# ------------------------------------------------------------------------------

env <- read.csv(ENV_RAW, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(env) < 2L) stop("Environmental table is empty or malformed: ", ENV_RAW)

# The first column of the uploaded environmental table stores the sample ID.
names(env)[1] <- "Sample_name"
env$Sample_name <- as.character(env$Sample_name)
env$plot <- as.character(env$plot)
env$block <- as.character(env$block)
env$year <- as.integer(env$year)

treat <- read.csv(TREAT_RAW, check.names = FALSE, stringsAsFactors = FALSE)
required_treat <- c("plot", "block", "Treatment", "wetting", "drying")
missing_treat <- setdiff(required_treat, names(treat))
if (length(missing_treat)) {
  stop("Treatment table is missing: ", paste(missing_treat, collapse = ", "))
}

treat$plot <- as.character(treat$plot)
treat$block <- as.character(treat$block)

# Add design variables required by the downstream split-plot scripts.
env$NS <- substr(env$plot, 1L, 1L)
env$Number <- sub("^[NS]", "", env$plot)
env$WholePlotID <- paste(env$block, env$Number, sep = "_")

env <- merge(
  env,
  treat[, c("plot", "Treatment", "wetting", "drying")],
  by = "plot",
  all.x = TRUE,
  sort = FALSE
)

# merge() can change row order; restore chronological/sample order.
row_order <- match(env$Sample_name, read.csv(ENV_RAW, check.names = FALSE,
                                              stringsAsFactors = FALSE)[[1]])
env <- env[order(row_order), , drop = FALSE]

if (anyNA(env$Treatment)) stop("At least one sample could not be matched to treatment_list10_20.csv")

# Consistent treatment labels used by the revised scripts.
env$Treatment_standard <- env$Treatment

env$WarmF <- factor(env$warm, levels = c(0, 1), labels = c("Unwarmed", "Warmed"))
env$PrecipF <- factor(env$precip, levels = c(1, 0.5, 2),
                      labels = c("Ambient", "Half", "Double"))
env$ClipF <- factor(env$clip, levels = c(0, 1), labels = c("Unclipped", "Clipped"))

# Keep a plain CSV with Sample_name as an ordinary column; this is the canonical
# metadata file used by all empirical scripts.
META_OUT <- file.path(PREP_DIR, "environment_metadata.csv")
write.csv(env, META_OUT, row.names = FALSE, quote = FALSE)

# 48-row split-plot metadata used by SEM and mixed-model scripts.
plot_meta <- treat
plot_meta$Number <- sub("^[NS]", "", plot_meta$plot)
plot_meta$NS <- substr(plot_meta$plot, 1L, 1L)
plot_meta$WholePlotID <- paste(plot_meta$block, plot_meta$Number, sep = "_")
PLOT_META_OUT <- file.path(PREP_DIR, "treatment_metadata_48plots.csv")
write.csv(plot_meta, PLOT_META_OUT, row.names = FALSE, quote = FALSE)

# Sample-level treatment table with Sample_name in the first column, suitable
# for iCAMP::match.name and icamp.bins.
icamp_treat <- env[, c("Sample_name", "Treatment", "Warm", "Precipitation", "Clip",
                       "warm", "wetting", "drying", "clip"), drop = FALSE]
write.csv(icamp_treat, file.path(PREP_DIR, "iCAMP_treatment_by_sample.csv"),
          row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# 2. Prepare plot-level relative-abundance matrices
# ------------------------------------------------------------------------------

group_files <- c(
  Archaea  = "resample7842_otutab_Archaea_2010to2020.csv",
  Bacteria = "resample26890_otutab_Bacteria_2010to2020.csv",
  Fungi    = "resample541_otutab_Fungi_2010to2020.csv",
  Protists = "resample1064_otutab_Protists_2010to2020.csv",
  Nematoda = "resample296_otutab_Nematoda_2010to2020.csv"
)

read_count_matrix <- function(file) {
  x <- read.csv(file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  x[] <- lapply(x, function(z) suppressWarnings(as.numeric(z)))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) stop("Non-numeric/NA entries in: ", file)
  if (any(x < 0)) stop("Negative abundance values in: ", file)
  x
}

for (group in names(group_files)) {
  message("Preparing ", group)
  infile <- file.path(DATA_DIR, group_files[[group]])
  if (!file.exists(infile)) stop("Missing community table: ", infile)

  counts <- read_count_matrix(infile)
  keep_samples <- intersect(colnames(counts), env$Sample_name)
  if (length(keep_samples) != nrow(env)) {
    missing <- setdiff(env$Sample_name, colnames(counts))
    stop(group, ": expected ", nrow(env), " sample columns but matched ",
         length(keep_samples), ". Missing examples: ",
         paste(head(missing, 10), collapse = ", "))
  }
  counts <- counts[, env$Sample_name, drop = FALSE]

  sample_totals <- colSums(counts)
  if (any(sample_totals <= 0)) stop(group, ": sample(s) with zero total abundance")
  rel <- sweep(counts, 2L, sample_totals, "/")

  full_out <- file.path(FULL_DIR, paste0(group, "_relative_abundance_2010to2020.csv"))
  write.csv(rel, full_out, quote = FALSE)

  gdir <- file.path(PLOT_DIR, group)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

  for (plot_id in unique(env$plot)) {
    md <- env[env$plot == plot_id, , drop = FALSE]
    md <- md[order(md$year), , drop = FALSE]
    if (!identical(md$year, 2010:2020)) {
      stop(group, " / ", plot_id, ": expected one sample for every year 2010-2020")
    }
    ids <- md$Sample_name
    plot_mat <- rel[, ids, drop = FALSE]
    write.csv(
      plot_mat,
      file.path(gdir, paste0(plot_id, "_UPARSE_otu_All_abun.csv")),
      quote = FALSE
    )
  }
}

# ------------------------------------------------------------------------------
# 3. Harmonize taxonomy files to a single TaxonID-first-column convention
# ------------------------------------------------------------------------------

tax_files <- c(
  Archaea  = "resample7842_otutab_Archaea_2010to2020_tax.csv",
  Bacteria = "resample26890_otutab_Bacteria_2010to2020_tax.csv",
  Fungi    = "resample541_otutab_Fungi_2010to2020_tax.csv",
  Protists = "resample1064_otutab_Protists_2010to2020_tax.csv",
  Nematoda = "resample296_otutab_Nematoda_2010to2020_tax_trophic.csv"
)

for (group in names(tax_files)) {
  infile <- file.path(DATA_DIR, tax_files[[group]])
  if (!file.exists(infile)) stop("Missing taxonomy file: ", infile)
  tx <- read.csv(infile, check.names = FALSE, stringsAsFactors = FALSE)

  if (group == "Nematoda" && "OTU" %in% names(tx)) {
    id <- as.character(tx$OTU)
  } else {
    id <- as.character(tx[[1]])
  }
  if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id)) {
    stop(group, ": taxonomy taxon IDs are missing or duplicated")
  }

  tx <- cbind(TaxonID = id, tx[, setdiff(seq_len(ncol(tx)), 1L), drop = FALSE])
  # Avoid duplicate TaxonID/OTU columns for Nematoda.
  if ("OTU" %in% names(tx)) tx$OTU <- NULL
  write.csv(tx, file.path(TAX_DIR, paste0(group, "_taxonomy.csv")),
            row.names = FALSE, quote = FALSE)
}

# ------------------------------------------------------------------------------
# 4. Input manifest
# ------------------------------------------------------------------------------

manifest <- data.frame(
  File = c(ENV_RAW, TREAT_RAW, file.path(DATA_DIR, unname(group_files)),
           file.path(DATA_DIR, unname(tax_files))),
  Exists = file.exists(c(ENV_RAW, TREAT_RAW, file.path(DATA_DIR, unname(group_files)),
                         file.path(DATA_DIR, unname(tax_files)))),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(PREP_DIR, "input_manifest.csv"), row.names = FALSE, quote = FALSE)

message("\nPreparation completed.")
message("Canonical metadata: ", META_OUT)
message("Plot-level community matrices: ", PLOT_DIR)
message("Harmonized taxonomy: ", TAX_DIR)
