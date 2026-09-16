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

## Test Taylor's power law (all ASVs/OTUs)
## tTTPL Top-10 process-intensity results

seed_list = c("0","1","2")
intensity_list = c(0, 0.25, 0.50, 0.75, 1.00)
process_list = c("immigration","growth","interaction","noise")
condition_list = c("warming","unwarming")

## Expected files for each process x condition:
## 10 Top-10 parameter sets x 5 intensities x 3 seeds = 150
EXPECTED_FILES_PER_PROCESS_CONDITION = 150
STOP_IF_INCOMPLETE = TRUE

##============================================================
## Parse the process-intensity filename.
##
## Example:
## Time_series_from_SLM_top10rank1_immigration_d0.0001_
## interaction_sd0.08_growth_g2_seed0_noise0.01_
## process_immigration_intensity0p25.csv
##============================================================

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

  m = regexec(pattern, basename(file_name))
  hit = regmatches(basename(file_name), m)[[1]]

  if(length(hit) != 9){
    stop(paste0("Cannot parse filename:\n", file_name))
  }

  intensity_value = as.numeric(gsub("p", ".", hit[9], fixed = TRUE))

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


##============================================================
## Calculate tTTPL beta for one resampled time-series file.
##
## This follows the previous min10 TPL code:
## - transpose the resampled table
## - use Y2010-Y2020 (first 11 time points)
## - calculate temporal mean and variance for each ASV/OTU
## - log10 transform mean and variance
## - fit lm(LogVar ~ LogMean)
##
## Variance > 0 is also required before log10 to prevent -Inf.
##============================================================

calculate_tTTPL_one <- function(data_name, summary_name){

  subcom = as.data.frame(
    t(read.csv(data_name, header = TRUE, row.names = 1, check.names = FALSE))
  )

  if(nrow(subcom) < 11){
    stop(
      paste0(
        "Fewer than 11 time points in:\n",
        data_name
      )
    )
  }

  ## Keep the same temporal window as the previous code:
  ## Y2010-Y2020; Y2021 is not used in the TPL regression.
  subcom = subcom[1:11, , drop = FALSE]

  all_summary = data.frame(
    Mean = double(),
    Variance = double()
  )

  for(j in 1:ncol(subcom)){

    abun = as.numeric(subcom[,j])

    mean_abun = mean(abun)
    var_abun = var(abun)

    abun_summary = data.frame(
      Mean = mean_abun,
      Variance = var_abun
    )

    rownames(abun_summary) = colnames(subcom)[j]

    all_summary = rbind(
      all_summary,
      abun_summary
    )
  }

  ## Mean > 0 follows the previous scripts.
  ## Variance > 0 is needed for a finite log10 variance.
  all_summary = all_summary[
    all_summary$Mean > 0 &
    all_summary$Variance > 0,
  ]

  if(nrow(all_summary) < 3){
    stop(
      paste0(
        "Too few valid ASVs/OTUs for TPL regression in:\n",
        data_name
      )
    )
  }

  all_summary["LogMean"] = log10(all_summary$Mean)
  all_summary["LogVar"] = log10(all_summary$Variance)

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


##============================================================
## Run one taxonomic group.
##
## Input:
## ./Resample/<Taxa>/tTTPL/process_intensity_top10/
##     <process>/<condition>/
##
## Output:
## ./Taylor/<Taxa>/tTTPL/process_intensity_top10/
##     <process>/<condition>/
##
## For every resampled time series, a mean-variance summary is saved.
## For every process x condition, one beta summary table is saved.
## One combined long-format beta table is also saved for each taxon.
##============================================================

run_tTTPL_process_intensity <- function(taxa){

  cat("\n============================================================\n")
  cat("Start tTTPL process-intensity:", taxa, "\n")
  cat("============================================================\n")

  taxon_beta_all = data.frame()

  for(process in process_list){

    for(condition in condition_list){

      input_dir = file.path(
        "./Resample",
        taxa,
        "tTTPL",
        "process_intensity_top10",
        process,
        condition
      )

      output_dir = file.path(
        "./Taylor",
        taxa,
        "tTTPL",
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

      if(!dir.exists(output_dir)){
        dir.create(
          output_dir,
          recursive = TRUE
        )
      }

      files = list.files(
        input_dir,
        pattern = "^Time_series_from_SLM_top10rank.*\\.csv$",
        full.names = TRUE
      )

      if(length(files) == 0){
        stop(
          paste0(
            "No process-intensity resample files found in:\n",
            input_dir
          )
        )
      }

      if(length(files) != EXPECTED_FILES_PER_PROCESS_CONDITION){

        msg = paste0(
          taxa, " | ", process, " | ", condition,
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

      ## Parse metadata first so files are processed in:
      ## Rank -> Intensity -> Seed order.
      metadata_all = do.call(
        rbind,
        lapply(files, parse_process_intensity_filename)
      )

      metadata_all$File = files

      metadata_all = metadata_all[
        order(
          metadata_all$Rank,
          metadata_all$Intensity,
          as.numeric(metadata_all$Seed)
        ),
      ]

      beta_summary_all = data.frame()

      for(k in 1:nrow(metadata_all)){

        meta = metadata_all[k,]

        ## Safety check: directory process should equal filename process.
        if(meta$Process != process){
          stop(
            paste0(
              "Process mismatch in filename:\n",
              meta$File
            )
          )
        }

        ## Safety check: only the expected intensity values.
        if(!any(abs(meta$Intensity - intensity_list) < 1e-12)){
          stop(
            paste0(
              "Unexpected process intensity ",
              meta$Intensity,
              " in:\n",
              meta$File
            )
          )
        }

        base_name = tools::file_path_sans_ext(
          basename(meta$File)
        )

        ## Keep the previous style of writing one mean/variance
        ## table for every TPL regression.
        summary_name = file.path(
          output_dir,
          paste0(
            sub("^Time_series_from_", "", base_name),
            "_all_species_mean_var_summary.csv"
          )
        )

        beta_result = calculate_tTTPL_one(
          data_name = meta$File,
          summary_name = summary_name
        )

        beta_summary = data.frame(
          Rank = meta$Rank,
          Immigration = meta$Immigration,
          Interaction = meta$Interaction,
          Growth = meta$Growth,
          Seed = meta$Seed,
          Noise = meta$Noise,
          Process = meta$Process,
          Intensity = meta$Intensity,
          Condition = condition,
          Intercept = beta_result$Intercept,
          Beta = beta_result$Beta,
          R2 = beta_result$R2,
          P = beta_result$P,
          stringsAsFactors = FALSE
        )

        beta_summary_all = rbind(
          beta_summary_all,
          beta_summary
        )

        cat(
          taxa,
          "|", process,
          "|", condition,
          "| rank", meta$Rank,
          "| intensity", meta$Intensity,
          "| seed", meta$Seed,
          "| beta", beta_result$Beta,
          "\n"
        )
      }

      beta_file = file.path(
        output_dir,
        paste0(
          "SLM_all_top10rank_immigration_interaction_growth_",
          "all_seed_noise_", process,
          "_process_intensity_all_species_Taylor_beta.csv"
        )
      )

      write.csv(
        beta_summary_all,
        beta_file,
        quote = FALSE,
        row.names = FALSE
      )

      taxon_beta_all = rbind(
        taxon_beta_all,
        beta_summary_all
      )
    }
  }

  ## Combined long-format result for downstream process-intensity
  ## curves and Top-10 robustness analyses.
  combined_dir = file.path(
    "./Taylor",
    taxa,
    "tTTPL",
    "process_intensity_top10"
  )

  if(!dir.exists(combined_dir)){
    dir.create(
      combined_dir,
      recursive = TRUE
    )
  }

  combined_file = file.path(
    combined_dir,
    "SLM_all_process_intensity_top10_all_species_Taylor_beta.csv"
  )

  write.csv(
    taxon_beta_all,
    combined_file,
    quote = FALSE,
    row.names = FALSE
  )

  cat("\nFinished:", taxa, "\n")
  cat("Combined beta table:", combined_file, "\n")
}


##============================================================
## Run all currently simulated tTTPL process-intensity taxa
##============================================================

taxa_list = c(
  "Archaea",
  "Bacteria",
  "Fungi_18S",
  "Protists",
  "Nematoda"
)

# taxa_list = c(
#   "Archaea"
# )

for(taxa in taxa_list){
  run_tTTPL_process_intensity(taxa)
}

cat("\nAll tTTPL process-intensity TPL tests finished.\n")
