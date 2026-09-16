rm(list=ls())
# ------------------------------------------------------------------------------
# PROJECT / sGLV ROOT
# ------------------------------------------------------------------------------
# Run from the repository root. Generated sGLV products are read from
# ./results/17_sGLV; no external input path outside ./data is required here.
PROJECT_DIR <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
SGLV_ROOT <- file.path(PROJECT_DIR, "results", "17_sGLV")
dir.create(SGLV_ROOT, recursive = TRUE, showWarnings = FALSE)
setwd(SGLV_ROOT)
PHYLO_ROOT <- file.path(PROJECT_DIR, "results", "pairwise_distance")

## ============================================================
## pTTPL Top-10 process-intensity Taylor-law analysis
##
## Input:
## ./Resample/<Taxa>/pTTPL/process_intensity_top10/
##     <process>/<condition>/
##
## MPD output:
## ./MPD/<Taxa>/pTTPL/process_intensity_top10/
##     <process>/<condition>/
##
## Taylor output:
## ./Taylor/<Taxa>/pTTPL/process_intensity_top10/
##     <process>/<condition>/
##
## The input filenames are produced by:
## 18_sGLV_resample.R
##
## For every process x condition:
## 10 Top-10 parameter sets x 5 intensities x 3 seeds = 150 files
## ============================================================

if (!requireNamespace("iCAMP", quietly = TRUE)) stop("Install R package iCAMP.")
if (!requireNamespace("bigmemory", quietly = TRUE)) stop("Install R package bigmemory.")

seed_list = c("0","1","2")
intensity_list = c(0, 0.25, 0.50, 0.75, 1.00)
process_list = c("immigration","growth","interaction","noise")
condition_list = c("warming","unwarming")

EXPECTED_FILES_PER_PROCESS_CONDITION = 150
STOP_IF_INCOMPLETE = TRUE


## ============================================================
## Calculate MPD_i,k of each OTU in each sample
##
## This follows the previous pTTPL script:
##   1. convert each sample to relative abundance
##   2. calculate abundance-weighted distance contribution
##   3. multiply by focal OTU relative abundance
##
## Final normalized taxon-specific contribution is obtained as:
##   MPD_i,k = each_otu_mpd / community_mpd
## ============================================================

mpdi.big <- function(comm,
                     pd.desc = "pd.desc",
                     pd.spname = NULL,
                     pd.wd = getwd(),
                     sp.limit = 10000,
                     spname.check = TRUE,
                     abundance.weighted = TRUE) {

  requireNamespace("bigmemory")

  if(spname.check){
    if(is.null(pd.spname)){
      pd.spname = read.csv(file.path(pd.wd, "pd.taxon.name.csv"), row.names = 1, stringsAsFactors = FALSE)[,1]
    }

    check.sp = iCAMP::match.name(
      name.check = pd.spname,
      cn.list = list(comm = comm)
    )

    comm = check.sp$comm
  }

  pdbig.id = match(colnames(comm), pd.spname)

  if(any(is.na(pdbig.id))){
    stop("Some community taxa could not be matched to the phylogenetic-distance matrix.")
  }

  if(!abundance.weighted){
    comm[comm > 0] = 1
  }

  rs = rowSums(comm)

  if(any(!is.finite(rs)) || any(rs <= 0)){
    stop("At least one sample has a non-positive or non-finite total abundance.")
  }

  comm = comm / rs
  comm = as.matrix(comm)

  pd = try(
    bigmemory::attach.big.matrix(
      dget(paste0(pd.wd, "/", pd.desc))
    ),
    silent = TRUE
  )

  if(inherits(pd, "try-error")){
    pd = bigmemory::attach.big.matrix(
      paste0(pd.wd, "/", pd.desc)
    )
  }

  ## All current microbial groups are far below sp.limit,
  ## but retain a chunked branch for completeness.
  if(nrow(pd) <= sp.limit){

    comd = comm %*% pd[pdbig.id, pdbig.id]

  } else {

    nsp = ncol(comm)
    N.int = max(1, floor(sp.limit * sp.limit / nsp))
    starts = seq(1, nsp, by = N.int)

    comd = matrix(
      0,
      nrow = nrow(comm),
      ncol = nsp
    )

    for(ii in seq_along(starts)){

      from = starts[ii]
      to = min(from + N.int - 1, nsp)

      comd[,from:to] =
        comm %*%
        pd[pdbig.id, pdbig.id[from:to], drop = FALSE]
    }
  }

  res = comd * comm
  colnames(res) = colnames(comm)

  res
}


## ============================================================
## Parse process-intensity filename
##
## Example:
## Time_series_from_SLM_top10rank1_immigration_d0.0001_
## interaction_sd0.4_growth_g3_seed0_noise0.05_
## process_immigration_intensity0p25.csv
## ============================================================

parse_process_intensity_filename <- function(file_name){

  pattern = paste0(
    "^Time_series_from_SLM_top10rank([0-9]+)",
    "_immigration_d([^_]+)",
    "_interaction_sd([^_]+)",
    "_growth_g([^_]+)",
    "_seed([^_]+)",
    "_noise([^_]+)",
    "_process_([^_]+)",
    "_intensity([^\\.]+)\\.csv$"
  )

  m = regexec(
    pattern,
    basename(file_name)
  )

  hit = regmatches(
    basename(file_name),
    m
  )[[1]]

  if(length(hit) != 9){
    stop(
      paste0(
        "Cannot parse filename:\n",
        file_name
      )
    )
  }

  intensity_value =
    as.numeric(
      gsub(
        "p",
        ".",
        hit[9],
        fixed = TRUE
      )
    )

  data.frame(
    Rank = as.integer(hit[2]),
    Immigration = as.numeric(hit[3]),
    Interaction = as.numeric(hit[4]),
    Growth = as.numeric(hit[5]),
    Seed = hit[6],
    Noise = as.numeric(hit[7]),
    Process = hit[8],
    Intensity = intensity_value,
    stringsAsFactors = FALSE
  )
}


## ============================================================
## Pairwise-distance files generated by 03_pTTPL_prepare_MPDik_by_plot.R
## ============================================================

taxa_config = list(
  Archaea = list(pdist_candidates = file.path(PHYLO_ROOT, "Archaea", "All_OTUs_2010to2020_pdist_out.rda")),
  Bacteria = list(pdist_candidates = file.path(PHYLO_ROOT, "Bacteria", "All_OTUs_2010to2020_pdist_out.rda")),
  Fungi_18S = list(pdist_candidates = file.path(PHYLO_ROOT, "Fungi", "All_OTUs_2010to2020_pdist_out.rda")),
  Protists = list(pdist_candidates = file.path(PHYLO_ROOT, "Protists", "All_OTUs_2010to2020_pdist_out.rda")),
  Nematoda = list(pdist_candidates = file.path(PHYLO_ROOT, "Nematoda", "All_OTUs_2010to2020_pdist_out.rda"))
)

resolve_pdist_file <- function(candidates, taxa){

  hit = candidates[file.exists(candidates)]

  if(length(hit) == 0){

    stop(
      paste0(
        "No phylogenetic-distance RDA file found for ",
        taxa,
        ". Checked:\n",
        paste(
          paste0("  ", candidates),
          collapse = "\n"
        )
      )
    )
  }

  hit[1]
}


load_pdist_out <- function(pdist_file){

  e = new.env(
    parent = emptyenv()
  )

  load(
    pdist_file,
    envir = e
  )

  if(!exists(
    "pdist.out",
    envir = e,
    inherits = FALSE
  )){
    stop(
      paste0(
        "Object 'pdist.out' was not found in:\n",
        pdist_file
      )
    )
  }

  obj = get(
    "pdist.out",
    envir = e,
    inherits = FALSE
  )

  ## Make copied/released pairwise-distance bundles portable.
  ## If the original absolute pd.wd no longer exists, use the directory that
  ## contains the RDA file, where pd.desc/pd.bin are expected to be released.
  desc_current = file.path(obj$pd.wd, obj$pd.file)
  if(!file.exists(desc_current)){
    desc_local = file.path(dirname(pdist_file), obj$pd.file)
    if(file.exists(desc_local)){
      obj$pd.wd = dirname(pdist_file)
    }else{
      stop(
        paste0(
          "Phylogenetic-distance descriptor was not found at either:\n  ",
          desc_current,
          "\n  ",
          desc_local
        )
      )
    }
  }

  obj
}


## ============================================================
## Calculate pTTPL beta for one resampled time-series file
##
## This follows the previous pTTPL workflow:
##   abundance time series
##       -> community MPD
##       -> taxon-specific MPD contribution
##       -> normalized MPD_i,k
##       -> temporal mean and variance across Y2010-Y2020
##       -> lm(LogVar ~ LogMean)
## ============================================================

calculate_pTTPL_one <- function(data_name,
                                mpd_name,
                                summary_name,
                                pdist.out){

  subcom = as.data.frame(
    t(
      read.csv(
        data_name,
        header = TRUE,
        row.names = 1,
        check.names = FALSE
      )
    )
  )

  if(nrow(subcom) < 11){
    stop(
      paste0(
        "Fewer than 11 time points in:\n",
        data_name
      )
    )
  }

  ## Same temporal window as the original tTTPL/pTTPL code:
  ## Y2010-Y2020; Y2021 is not used in the Taylor-law fit.
  subcom = subcom[
    1:11,
    ,
    drop = FALSE
  ]

  ## Community-level abundance-weighted MPD
  subcom_mpd = iCAMP::mpd.big(
    subcom,
    pd.desc = pdist.out$pd.file,
    pd.spname = pdist.out$tip.label,
    pd.wd = pdist.out$pd.wd,
    spname.check = TRUE
  )

  ## Taxon-specific unnormalized MPD contribution
  each_otu_mpd = mpdi.big(
    subcom,
    pd.desc = pdist.out$pd.file,
    pd.spname = pdist.out$tip.label,
    pd.wd = pdist.out$pd.wd,
    spname.check = TRUE
  )

  ## Same normalization as the previous pTTPL script
  otu_MPDik = each_otu_mpd / subcom_mpd

  if(any(!is.finite(otu_MPDik))){
    stop(
      paste0(
        "Non-finite MPD_i,k values produced for:\n",
        data_name
      )
    )
  }

  write.csv(
    otu_MPDik,
    mpd_name,
    quote = FALSE
  )

  all_summary = data.frame(
    Mean = double(),
    Variance = double()
  )

  for(j in 1:ncol(otu_MPDik)){

    abun = as.numeric(
      otu_MPDik[,j]
    )

    mean_abun = mean(abun)
    var_abun = var(abun)

    abun_summary = data.frame(
      Mean = mean_abun,
      Variance = var_abun
    )

    rownames(abun_summary) =
      colnames(otu_MPDik)[j]

    all_summary = rbind(
      all_summary,
      abun_summary
    )
  }

  ## Mean > 0 follows the previous pTTPL script.
  ## Variance > 0 additionally prevents log10(0) = -Inf.
  all_summary = all_summary[
    all_summary$Mean > 0 &
    all_summary$Variance > 0,
  ]

  if(nrow(all_summary) < 3){
    stop(
      paste0(
        "Too few valid taxa for the pTTPL regression in:\n",
        data_name
      )
    )
  }

  all_summary["LogMean"] =
    log10(
      all_summary$Mean
    )

  all_summary["LogVar"] =
    log10(
      all_summary$Variance
    )

  write.csv(
    all_summary,
    summary_name,
    quote = FALSE
  )

  model.lm = lm(
    LogVar ~ LogMean,
    data = all_summary
  )

  ss = summary(model.lm)

  ## Preserve full numerical precision in downstream sensitivity calculations.
  data.frame(
    Intercept = unname(coef(model.lm)[1]),
    Beta = unname(coef(model.lm)[2]),
    R2 = unname(ss$r.squared),
    P = unname(ss$coefficients[2,4])
  )
}


## ============================================================
## Run one taxonomic group
## ============================================================

run_pTTPL_process_intensity <- function(taxa,
                                        config){

  cat(
    "\n============================================================\n"
  )
  cat(
    "Start pTTPL process-intensity:",
    taxa,
    "\n"
  )
  cat(
    "============================================================\n"
  )

  pdist_file =
    resolve_pdist_file(
      config$pdist_candidates,
      taxa
    )

  cat(
    "Using phylogenetic-distance file:\n",
    pdist_file,
    "\n"
  )

  ## Load once per taxon, rather than once for every simulated file.
  pdist.out =
    load_pdist_out(
      pdist_file
    )

  taxon_beta_all =
    data.frame()

  for(process in process_list){

    for(condition in condition_list){

      input_dir = file.path(
        "./Resample",
        taxa,
        "pTTPL",
        "process_intensity_top10",
        process,
        condition
      )

      mpd_dir = file.path(
        "./MPD",
        taxa,
        "pTTPL",
        "process_intensity_top10",
        process,
        condition
      )

      output_dir = file.path(
        "./Taylor",
        taxa,
        "pTTPL",
        "process_intensity_top10",
        process,
        condition
      )

      if(!dir.exists(input_dir)){
        stop(
          paste0(
            "Input directory does not exist:\n",
            input_dir
          )
        )
      }

      if(!dir.exists(mpd_dir)){
        dir.create(
          mpd_dir,
          recursive = TRUE
        )
      }

      if(!dir.exists(output_dir)){
        dir.create(
          output_dir,
          recursive = TRUE
        )
      }

      files = list.files(
        input_dir,
        pattern =
          "^Time_series_from_SLM_top10rank.*\\.csv$",
        full.names = TRUE
      )

      if(length(files) == 0){
        stop(
          paste0(
            "No pTTPL process-intensity resample files found in:\n",
            input_dir
          )
        )
      }

      if(
        length(files) !=
        EXPECTED_FILES_PER_PROCESS_CONDITION
      ){

        msg = paste0(
          taxa,
          " | ",
          process,
          " | ",
          condition,
          ": expected ",
          EXPECTED_FILES_PER_PROCESS_CONDITION,
          " files but found ",
          length(files),
          "."
        )

        if(STOP_IF_INCOMPLETE){
          stop(msg)
        }else{
          warning(msg)
        }
      }

      ## Parse metadata so files are processed in:
      ## Rank -> Intensity -> Seed order.
      metadata_all =
        do.call(
          rbind,
          lapply(
            files,
            parse_process_intensity_filename
          )
        )

      metadata_all$File =
        files

      metadata_all =
        metadata_all[
          order(
            metadata_all$Rank,
            metadata_all$Intensity,
            as.numeric(
              metadata_all$Seed
            )
          ),
        ]

      beta_summary_all =
        data.frame()

      for(k in 1:nrow(metadata_all)){

        meta =
          metadata_all[k,]

        if(meta$Process != process){
          stop(
            paste0(
              "Process mismatch in filename:\n",
              meta$File
            )
          )
        }

        if(
          !any(
            abs(
              meta$Intensity -
              intensity_list
            ) < 1e-12
          )
        ){
          stop(
            paste0(
              "Unexpected process intensity ",
              meta$Intensity,
              " in:\n",
              meta$File
            )
          )
        }

        base_name =
          tools::file_path_sans_ext(
            basename(
              meta$File
            )
          )

        ## Example:
        ## MPD_weighted_SLM_top10rank1_..._intensity0p25.csv
        mpd_name =
          file.path(
            mpd_dir,
            paste0(
              sub(
                "^Time_series_from_",
                "MPD_weighted_",
                base_name
              ),
              ".csv"
            )
          )

        summary_name =
          file.path(
            output_dir,
            paste0(
              sub(
                "^Time_series_from_",
                "MPD_weighted_",
                base_name
              ),
              "_all_species_mean_var_summary.csv"
            )
          )

        beta_result =
          calculate_pTTPL_one(
            data_name = meta$File,
            mpd_name = mpd_name,
            summary_name = summary_name,
            pdist.out = pdist.out
          )

        beta_summary =
          data.frame(
            Rank = meta$Rank,
            Immigration =
              meta$Immigration,
            Interaction =
              meta$Interaction,
            Growth =
              meta$Growth,
            Seed =
              meta$Seed,
            Noise =
              meta$Noise,
            Process =
              meta$Process,
            Intensity =
              meta$Intensity,
            Condition =
              condition,
            Intercept =
              beta_result$Intercept,
            Beta =
              beta_result$Beta,
            R2 =
              beta_result$R2,
            P =
              beta_result$P,
            stringsAsFactors = FALSE
          )

        beta_summary_all =
          rbind(
            beta_summary_all,
            beta_summary
          )

        cat(
          taxa,
          "|",
          process,
          "|",
          condition,
          "| rank",
          meta$Rank,
          "| intensity",
          meta$Intensity,
          "| seed",
          meta$Seed,
          "| beta_p",
          beta_result$Beta,
          "\n"
        )
      }

      beta_file =
        file.path(
          output_dir,
          paste0(
            "MPD_weighted_SLM_all_top10rank_",
            "immigration_interaction_growth_",
            "all_seed_noise_",
            process,
            "_process_intensity_all_species_",
            "Taylor_beta.csv"
          )
        )

      write.csv(
        beta_summary_all,
        beta_file,
        quote = FALSE,
        row.names = FALSE
      )

      taxon_beta_all =
        rbind(
          taxon_beta_all,
          beta_summary_all
        )
    }
  }

  ## Combined long-format result for downstream
  ## process-intensity curves and Top-10 robustness analysis.
  combined_dir =
    file.path(
      "./Taylor",
      taxa,
      "pTTPL",
        "process_intensity_top10"
    )

  if(!dir.exists(combined_dir)){
    dir.create(
      combined_dir,
      recursive = TRUE
    )
  }

  combined_file =
    file.path(
      combined_dir,
      paste0(
        "MPD_weighted_SLM_all_process_intensity_",
        "top10_all_species_Taylor_beta.csv"
      )
    )

  write.csv(
    taxon_beta_all,
    combined_file,
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    "\nFinished:",
    taxa,
    "\n"
  )

  cat(
    "Combined beta_p table:",
    combined_file,
    "\n"
  )
}


## ============================================================
## Run all five taxa
## ============================================================

taxa_list = c(
  "Archaea",
  "Bacteria",
  "Fungi_18S",
  "Protists",
  "Nematoda"
)

## For a test run on Archaea only, use:
# taxa_list = c("Archaea")

for(taxa in taxa_list){

  run_pTTPL_process_intensity(
    taxa = taxa,
    config = taxa_config[[taxa]]
  )
}

cat(
  "\nAll pTTPL Top-10 process-intensity Taylor-law tests finished.\n"
)
