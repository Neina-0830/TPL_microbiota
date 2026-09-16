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

library(dplyr)
library(tidyr)
library(ggplot2)

## Use explicit namespaces for dplyr/tidyr verbs below to avoid
## conflicts such as MASS::select masking dplyr::select in an existing R session.

## ============================================================
## Run process-intensity importance analysis for:
##   5 taxa x (tTTPL + pTTPL)
##
## tTTPL input:
## ./Taylor/<Taxa>/tTTPL/process_intensity_top10/
##   SLM_all_process_intensity_top10_all_species_Taylor_beta.csv
##
## pTTPL input:
## ./Taylor/<Taxa>/pTTPL/process_intensity_top10/
##   MPD_weighted_SLM_all_process_intensity_top10_all_species_Taylor_beta.csv
##
## Per-taxon outputs are stored in:
##   <combined_dir>/evaluation/
##
## Cross-taxa combined outputs are stored in:
## ./Taylor/process_intensity_top10_all_taxa_evaluation/
## ============================================================

taxa_list <- c("Bacteria","Archaea","Fungi_18S","Protists","Nematoda")
tpl_list <- c("tTTPL","pTTPL")
process_order <- c("immigration","growth","interaction","noise")
condition_order <- c("warming","unwarming")
intensity_order <- c(0,0.25,0.50,0.75,1.00)

EXPECTED_ROWS <- 1200
STOP_IF_INCOMPLETE <- TRUE

## BUGFIX NOTE:
## Relative importance is now normalized with explicit group-wise
## denominators (AUC_total and Endpoint_total). This avoids the prior
## scalar-ifelse recycling bug that assigned the same RI to all processes.


## ---------- path helpers ----------

get_combined_dir <- function(taxa, tpl){
  if(tpl == "tTTPL"){
    file.path("./Taylor", taxa, "tTTPL", "process_intensity_top10")
  } else {
    file.path("./Taylor", taxa, "pTTPL", "process_intensity_top10")
  }
}

get_input_file <- function(taxa, tpl){
  d <- get_combined_dir(taxa, tpl)
  if(tpl == "tTTPL"){
    file.path(d, "SLM_all_process_intensity_top10_all_species_Taylor_beta.csv")
  } else {
    file.path(d, "MPD_weighted_SLM_all_process_intensity_top10_all_species_Taylor_beta.csv")
  }
}


## ---------- helper functions ----------

trapz_auc <- function(x, y){
  o <- order(x)
  x <- x[o]
  y <- y[o]
  sum(diff(x) * (head(y,-1) + tail(y,-1))/2)
}

safe_spearman <- function(x, y){
  if(length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(cor(x, y, method="spearman", use="complete.obs"))
}

q25 <- function(x) quantile(x, 0.25, na.rm=TRUE, names=FALSE)
q75 <- function(x) quantile(x, 0.75, na.rm=TRUE, names=FALSE)


## ============================================================
## Evaluate one Taxa x TPL combination
## ============================================================

evaluate_one <- function(taxa, tpl){
  
  input_file <- get_input_file(taxa, tpl)
  combined_dir <- get_combined_dir(taxa, tpl)
  output_dir <- file.path(combined_dir, "evaluation")
  
  cat("\n--------------------------------------------------\n")
  cat("Running:", taxa, "|", tpl, "\n")
  cat("Input:", input_file, "\n")
  
  if(!file.exists(input_file)){
    stop(paste0("Input file not found:\n", input_file))
  }
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=TRUE)
  }
  
  dat <- read.csv(
    input_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
  ) %>%
    dplyr::mutate(
      Rank=as.integer(Rank),
      Seed=as.character(Seed),
      Intensity=as.numeric(Intensity),
      Beta=as.numeric(Beta),
      Process=as.character(Process),
      Condition=as.character(Condition)
    )
  
  required <- c(
    "Rank","Immigration","Interaction","Growth","Seed","Noise",
    "Process","Intensity","Condition","Intercept","Beta","R2","P"
  )
  
  miss <- setdiff(required, colnames(dat))
  if(length(miss)>0){
    stop(paste("Missing columns:", paste(miss, collapse=", ")))
  }
  
  if(nrow(dat) != EXPECTED_ROWS){
    msg <- paste0(taxa," ",tpl,": expected ",EXPECTED_ROWS,
                  " rows, found ",nrow(dat))
    if(STOP_IF_INCOMPLETE) stop(msg) else warning(msg)
  }
  
  if(!setequal(unique(dat$Process), process_order)){
    stop(paste(taxa,tpl,": unexpected Process values"))
  }
  
  if(!setequal(unique(dat$Condition), condition_order)){
    stop(paste(taxa,tpl,": unexpected Condition values"))
  }
  
  if(!all(sort(unique(dat$Intensity)) == intensity_order)){
    stop(paste(taxa,tpl,": unexpected intensity grid"))
  }
  
  dup_check <- dat %>%
    dplyr::count(Condition,Process,Rank,Seed,Intensity) %>%
    dplyr::filter(n != 1)
  
  if(nrow(dup_check)>0){
    stop(paste(taxa,tpl,": duplicated/incomplete grid"))
  }
  
  
  ## 1. z=1 QC
  qc_z1 <- dat %>%
    dplyr::filter(Intensity==1) %>%
    dplyr::select(Condition,Rank,Seed,Process,Beta) %>%
    tidyr::pivot_wider(names_from=Process, values_from=Beta,
                       names_prefix="Beta_z1_") %>%
    rowwise() %>%
    dplyr::mutate(
      Beta_z1_mean=mean(c_across(starts_with("Beta_z1_")),na.rm=TRUE),
      Beta_z1_SD=sd(c_across(starts_with("Beta_z1_")),na.rm=TRUE),
      Beta_z1_range=
        max(c_across(starts_with("Beta_z1_")),na.rm=TRUE) -
        min(c_across(starts_with("Beta_z1_")),na.rm=TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(qc_z1,
            file.path(output_dir,"01_QC_z1_full_model_across_processes.csv"),
            row.names=FALSE)
  
  
  ## 2. seed-level curves
  seed_curve <- dat %>%
    dplyr::group_by(Condition,Rank,Seed,Process) %>%
    dplyr::arrange(Intensity,.by_group=TRUE) %>%
    dplyr::mutate(
      BetaFull_z1=Beta[Intensity==1][1],
      Effect_to_full=BetaFull_z1-Beta,
      AbsEffect_to_full=abs(Effect_to_full)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(seed_curve,
            file.path(output_dir,"02_beta_intensity_seed_level.csv"),
            row.names=FALSE)
  
  
  ## 3. seed-level AUC / endpoint / monotonicity
  seed_metric <- seed_curve %>%
    dplyr::group_by(Condition,Rank,Seed,Process) %>%
    dplyr::arrange(Intensity,.by_group=TRUE) %>%
    dplyr::summarise(
      BetaFull_z1=Beta[Intensity==1][1],
      AUC_abs=trapz_auc(Intensity,abs(Effect_to_full)),
      AUC_signed=trapz_auc(Intensity,Effect_to_full),
      Endpoint_abs=abs(Beta[Intensity==1][1]-Beta[Intensity==0][1]),
      Endpoint_signed=Beta[Intensity==1][1]-Beta[Intensity==0][1],
      Spearman_rho=safe_spearman(Intensity,Beta),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(seed_metric,
            file.path(output_dir,"03_seed_level_sensitivity_metrics.csv"),
            row.names=FALSE)
  
  
  ## 4. average 3 seeds within rank: curve table
  rank_curve <- seed_curve %>%
    dplyr::group_by(Condition,Rank,Process,Intensity) %>%
    dplyr::summarise(
      Beta_mean_seed=mean(Beta,na.rm=TRUE),
      Beta_SD_seed=sd(Beta,na.rm=TRUE),
      Effect_to_full_mean_seed=mean(Effect_to_full,na.rm=TRUE),
      n_seed=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(rank_curve,
            file.path(output_dir,"04_beta_intensity_rank_mean.csv"),
            row.names=FALSE)
  
  
  ## 5. plot-ready Top10 beta curve summary
  curve_summary <- rank_curve %>%
    dplyr::group_by(Condition,Process,Intensity) %>%
    dplyr::summarise(
      Beta_mean_Top10=mean(Beta_mean_seed,na.rm=TRUE),
      Beta_median_Top10=median(Beta_mean_seed,na.rm=TRUE),
      Beta_SD_Top10=sd(Beta_mean_seed,na.rm=TRUE),
      Beta_Q25_Top10=q25(Beta_mean_seed),
      Beta_Q75_Top10=q75(Beta_mean_seed),
      Effect_to_full_mean_Top10=
        mean(Effect_to_full_mean_seed,na.rm=TRUE),
      n_rank=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(curve_summary,
            file.path(output_dir,"05_beta_intensity_Top10_summary_for_plot.csv"),
            row.names=FALSE)
  
  
  ## 6. rank-level process importance
  ## first average 3 seeds; then normalize four processes within Rank
  importance_rank <- seed_metric %>%
    dplyr::group_by(Condition,Rank,Process) %>%
    dplyr::summarise(
      AUC_abs=mean(AUC_abs,na.rm=TRUE),
      AUC_signed=mean(AUC_signed,na.rm=TRUE),
      Endpoint_abs=mean(Endpoint_abs,na.rm=TRUE),
      Endpoint_signed=mean(Endpoint_signed,na.rm=TRUE),
      Spearman_rho=mean(Spearman_rho,na.rm=TRUE),
      BetaFull_z1_mean_seed=mean(BetaFull_z1,na.rm=TRUE),
      n_seed=n(),
      .groups="drop"
    ) %>%
    dplyr::group_by(Condition,Rank) %>%
    dplyr::mutate(
      AUC_total =
        sum(AUC_abs, na.rm = TRUE),
      
      Endpoint_total =
        sum(Endpoint_abs, na.rm = TRUE),
      
      RelativeImportance_AUC =
        dplyr::if_else(
          AUC_total > 0,
          AUC_abs / AUC_total,
          NA_real_
        ),
      
      RelativeImportance_Endpoint =
        dplyr::if_else(
          Endpoint_total > 0,
          Endpoint_abs / Endpoint_total,
          NA_real_
        ),
      
      Top1_AUC =
        as.integer(
          RelativeImportance_AUC ==
            max(
              RelativeImportance_AUC,
              na.rm = TRUE
            )
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  ## QC: within every Condition x Rank, the four process RIs
  ## must sum to 1 (unless all AUCs are zero).
  ri_qc <- importance_rank %>%
    dplyr::group_by(Condition, Rank) %>%
    dplyr::summarise(
      AUC_total_check = dplyr::first(AUC_total),
      RI_AUC_sum = sum(RelativeImportance_AUC, na.rm = TRUE),
      RI_Endpoint_sum = sum(RelativeImportance_Endpoint, na.rm = TRUE),
      Top1_count = sum(Top1_AUC, na.rm = TRUE),
      .groups = "drop"
    )
  
  bad_ri <- ri_qc %>%
    dplyr::filter(
      AUC_total_check > 0 &
        abs(RI_AUC_sum - 1) > 1e-10
    )
  
  if(nrow(bad_ri) > 0){
    stop(
      paste0(
        taxa, " ", tpl,
        ": relative AUC importance does not sum to 1 within some Condition x Rank groups."
      )
    )
  }
  
  write.csv(
    ri_qc,
    file.path(
      output_dir,
      "06b_process_importance_normalization_QC.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(importance_rank,
            file.path(output_dir,"06_process_importance_rank_level.csv"),
            row.names=FALSE)
  
  
  ## 7. Top10 ensemble summary
  importance_summary <- importance_rank %>%
    dplyr::group_by(Condition,Process) %>%
    dplyr::summarise(
      Mean_RelativeImportance_AUC=
        mean(RelativeImportance_AUC,na.rm=TRUE),
      Median_RelativeImportance_AUC=
        median(RelativeImportance_AUC,na.rm=TRUE),
      SD_RelativeImportance_AUC=
        sd(RelativeImportance_AUC,na.rm=TRUE),
      Q25_RelativeImportance_AUC=q25(RelativeImportance_AUC),
      Q75_RelativeImportance_AUC=q75(RelativeImportance_AUC),
      Mean_AUC_abs=mean(AUC_abs,na.rm=TRUE),
      Mean_AUC_signed=mean(AUC_signed,na.rm=TRUE),
      Mean_RelativeImportance_Endpoint=
        mean(RelativeImportance_Endpoint,na.rm=TRUE),
      Mean_Endpoint_signed=mean(Endpoint_signed,na.rm=TRUE),
      Mean_Spearman_rho=mean(Spearman_rho,na.rm=TRUE),
      Top1_Count_out_of_10=sum(Top1_AUC,na.rm=TRUE),
      Positive_AUC_signed_Count=sum(AUC_signed>0,na.rm=TRUE),
      Negative_AUC_signed_Count=sum(AUC_signed<0,na.rm=TRUE),
      n_rank=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(importance_summary,
            file.path(output_dir,"07_process_importance_Top10_summary_for_plot.csv"),
            row.names=FALSE)
  
  
  ## 8. rank stability
  rank_stability <- importance_rank %>%
    dplyr::group_by(Condition,Rank) %>%
    dplyr::arrange(desc(RelativeImportance_AUC),.by_group=TRUE) %>%
    dplyr::mutate(Importance_Rank=row_number()) %>%
    dplyr::ungroup()
  
  write.csv(rank_stability,
            file.path(output_dir,"08_process_importance_rank_stability.csv"),
            row.names=FALSE)
  
  
  ## 9. warming-unwarming beta contrast
  warming_seed <- dat %>%
    dplyr::select(Rank,Seed,Process,Intensity,Condition,Beta) %>%
    tidyr::pivot_wider(names_from=Condition,values_from=Beta) %>%
    dplyr::mutate(
      DeltaBeta_WminusU=warming-unwarming
    )
  
  warming_rank_curve <- warming_seed %>%
    dplyr::group_by(Rank,Process,Intensity) %>%
    dplyr::summarise(
      DeltaBeta_mean_seed=mean(DeltaBeta_WminusU,na.rm=TRUE),
      DeltaBeta_SD_seed=sd(DeltaBeta_WminusU,na.rm=TRUE),
      n_seed=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(warming_rank_curve,
            file.path(output_dir,"09_warming_contrast_intensity_rank_mean.csv"),
            row.names=FALSE)
  
  warming_curve_summary <- warming_rank_curve %>%
    dplyr::group_by(Process,Intensity) %>%
    dplyr::summarise(
      DeltaBeta_mean_Top10=mean(DeltaBeta_mean_seed,na.rm=TRUE),
      DeltaBeta_median_Top10=median(DeltaBeta_mean_seed,na.rm=TRUE),
      DeltaBeta_SD_Top10=sd(DeltaBeta_mean_seed,na.rm=TRUE),
      DeltaBeta_Q25_Top10=q25(DeltaBeta_mean_seed),
      DeltaBeta_Q75_Top10=q75(DeltaBeta_mean_seed),
      n_rank=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(warming_curve_summary,
            file.path(output_dir,"10_warming_contrast_intensity_Top10_summary_for_plot.csv"),
            row.names=FALSE)
  
  
  ## 10. warming-contrast AUC importance
  warming_seed_metric <- warming_seed %>%
    dplyr::group_by(Rank,Seed,Process) %>%
    dplyr::arrange(Intensity,.by_group=TRUE) %>%
    dplyr::mutate(
      FullDeltaBeta_z1=DeltaBeta_WminusU[Intensity==1][1],
      Effect_to_full=FullDeltaBeta_z1-DeltaBeta_WminusU
    ) %>%
    dplyr::summarise(
      FullDeltaBeta_z1=first(FullDeltaBeta_z1),
      AUC_abs=trapz_auc(Intensity,abs(Effect_to_full)),
      AUC_signed=trapz_auc(Intensity,Effect_to_full),
      Endpoint_abs=abs(
        DeltaBeta_WminusU[Intensity==1][1] -
          DeltaBeta_WminusU[Intensity==0][1]
      ),
      Endpoint_signed=
        DeltaBeta_WminusU[Intensity==1][1] -
        DeltaBeta_WminusU[Intensity==0][1],
      Spearman_rho=safe_spearman(Intensity,DeltaBeta_WminusU),
      .groups="drop"
    )
  
  warming_importance_rank <- warming_seed_metric %>%
    dplyr::group_by(Rank,Process) %>%
    dplyr::summarise(
      AUC_abs=mean(AUC_abs,na.rm=TRUE),
      AUC_signed=mean(AUC_signed,na.rm=TRUE),
      Endpoint_abs=mean(Endpoint_abs,na.rm=TRUE),
      Endpoint_signed=mean(Endpoint_signed,na.rm=TRUE),
      Spearman_rho=mean(Spearman_rho,na.rm=TRUE),
      FullDeltaBeta_z1_mean_seed=mean(FullDeltaBeta_z1,na.rm=TRUE),
      n_seed=n(),
      .groups="drop"
    ) %>%
    dplyr::group_by(Rank) %>%
    dplyr::mutate(
      AUC_total =
        sum(AUC_abs, na.rm = TRUE),
      
      Endpoint_total =
        sum(Endpoint_abs, na.rm = TRUE),
      
      RelativeImportance_AUC =
        dplyr::if_else(
          AUC_total > 0,
          AUC_abs / AUC_total,
          NA_real_
        ),
      
      RelativeImportance_Endpoint =
        dplyr::if_else(
          Endpoint_total > 0,
          Endpoint_abs / Endpoint_total,
          NA_real_
        ),
      
      Top1_AUC =
        as.integer(
          RelativeImportance_AUC ==
            max(
              RelativeImportance_AUC,
              na.rm = TRUE
            )
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  ## QC: within every Rank, the four warming-contrast process RIs
  ## must sum to 1 (unless all AUCs are zero).
  warming_ri_qc <- warming_importance_rank %>%
    dplyr::group_by(Rank) %>%
    dplyr::summarise(
      AUC_total_check = dplyr::first(AUC_total),
      RI_AUC_sum = sum(RelativeImportance_AUC, na.rm = TRUE),
      RI_Endpoint_sum = sum(RelativeImportance_Endpoint, na.rm = TRUE),
      Top1_count = sum(Top1_AUC, na.rm = TRUE),
      .groups = "drop"
    )
  
  bad_warming_ri <- warming_ri_qc %>%
    dplyr::filter(
      AUC_total_check > 0 &
        abs(RI_AUC_sum - 1) > 1e-10
    )
  
  if(nrow(bad_warming_ri) > 0){
    stop(
      paste0(
        taxa, " ", tpl,
        ": warming-contrast relative AUC importance does not sum to 1 within some Rank groups."
      )
    )
  }
  
  write.csv(
    warming_ri_qc,
    file.path(
      output_dir,
      "11b_warming_contrast_importance_normalization_QC.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(warming_importance_rank,
            file.path(output_dir,
                      "11_warming_contrast_process_importance_rank_level.csv"),
            row.names=FALSE)
  
  warming_importance_summary <- warming_importance_rank %>%
    dplyr::group_by(Process) %>%
    dplyr::summarise(
      Mean_RelativeImportance_AUC=
        mean(RelativeImportance_AUC,na.rm=TRUE),
      Median_RelativeImportance_AUC=
        median(RelativeImportance_AUC,na.rm=TRUE),
      SD_RelativeImportance_AUC=
        sd(RelativeImportance_AUC,na.rm=TRUE),
      Q25_RelativeImportance_AUC=q25(RelativeImportance_AUC),
      Q75_RelativeImportance_AUC=q75(RelativeImportance_AUC),
      Mean_AUC_abs=mean(AUC_abs,na.rm=TRUE),
      Mean_AUC_signed=mean(AUC_signed,na.rm=TRUE),
      Mean_RelativeImportance_Endpoint=
        mean(RelativeImportance_Endpoint,na.rm=TRUE),
      Mean_Endpoint_signed=mean(Endpoint_signed,na.rm=TRUE),
      Mean_Spearman_rho=mean(Spearman_rho,na.rm=TRUE),
      Top1_Count_out_of_10=sum(Top1_AUC,na.rm=TRUE),
      n_rank=n(),
      .groups="drop"
    ) %>%
    dplyr::mutate(Taxa=taxa,TPL=tpl)
  
  write.csv(warming_importance_summary,
            file.path(output_dir,
                      "12_warming_contrast_process_importance_Top10_summary_for_plot.csv"),
            row.names=FALSE)
  
  cat("Finished:", taxa, "|", tpl, "\n")
  
  list(
    qc_z1=qc_z1,
    rank_curve=rank_curve,
    curve_summary=curve_summary,
    importance_rank=importance_rank,
    importance_summary=importance_summary,
    rank_stability=rank_stability,
    warming_rank_curve=warming_rank_curve,
    warming_curve_summary=warming_curve_summary,
    warming_importance_rank=warming_importance_rank,
    warming_importance_summary=warming_importance_summary
  )
}


## ============================================================
## Run 5 taxa x 2 TPL types
## ============================================================

results_all <- list()

for(tpl in tpl_list){
  for(taxa in taxa_list){
    key <- paste(tpl,taxa,sep="__")
    results_all[[key]] <- evaluate_one(taxa,tpl)
  }
}


## ============================================================
## Combine and save cross-taxa results
## ============================================================

all_out <- "./Taylor/process_intensity_top10_all_taxa_evaluation"

if(!dir.exists(all_out)){
  dir.create(all_out,recursive=TRUE)
}

combine_component <- function(name){
  dplyr::bind_rows(lapply(results_all,function(x)x[[name]]))
}

all_qc <- combine_component("qc_z1")
all_rank_curve <- combine_component("rank_curve")
all_curve_summary <- combine_component("curve_summary")
all_imp_rank <- combine_component("importance_rank")
all_imp_summary <- combine_component("importance_summary")
all_rank_stability <- combine_component("rank_stability")
all_warm_rank_curve <- combine_component("warming_rank_curve")
all_warm_curve_summary <- combine_component("warming_curve_summary")
all_warm_imp_rank <- combine_component("warming_importance_rank")
all_warm_imp_summary <- combine_component("warming_importance_summary")

write.csv(all_qc,
          file.path(all_out,"01_ALL_TAXA_QC_z1.csv"),
          row.names=FALSE)

write.csv(all_rank_curve,
          file.path(all_out,"02_ALL_TAXA_beta_intensity_rank_mean.csv"),
          row.names=FALSE)

write.csv(all_curve_summary,
          file.path(all_out,
                    "03_ALL_TAXA_beta_intensity_Top10_summary_for_plot.csv"),
          row.names=FALSE)

write.csv(all_imp_rank,
          file.path(all_out,
                    "04_ALL_TAXA_process_importance_rank_level.csv"),
          row.names=FALSE)

write.csv(all_imp_summary,
          file.path(all_out,
                    "05_ALL_TAXA_process_importance_Top10_summary_for_plot.csv"),
          row.names=FALSE)

write.csv(all_rank_stability,
          file.path(all_out,
                    "06_ALL_TAXA_process_importance_rank_stability.csv"),
          row.names=FALSE)

write.csv(all_warm_rank_curve,
          file.path(all_out,
                    "07_ALL_TAXA_warming_contrast_intensity_rank_mean.csv"),
          row.names=FALSE)

write.csv(all_warm_curve_summary,
          file.path(all_out,
                    "08_ALL_TAXA_warming_contrast_intensity_Top10_summary_for_plot.csv"),
          row.names=FALSE)

write.csv(all_warm_imp_rank,
          file.path(all_out,
                    "09_ALL_TAXA_warming_contrast_process_importance_rank_level.csv"),
          row.names=FALSE)

write.csv(all_warm_imp_summary,
          file.path(all_out,
                    "10_ALL_TAXA_warming_contrast_process_importance_Top10_summary_for_plot.csv"),
          row.names=FALSE)


## ============================================================
## Concise ranking tables
## ============================================================

condition_ranking <- all_imp_summary %>%
  dplyr::group_by(TPL,Taxa,Condition) %>%
  dplyr::arrange(desc(Mean_RelativeImportance_AUC),.by_group=TRUE) %>%
  dplyr::mutate(Mean_Importance_Rank=row_number()) %>%
  dplyr::ungroup()

write.csv(condition_ranking,
          file.path(all_out,
                    "11_ALL_TAXA_condition_specific_process_ranking.csv"),
          row.names=FALSE)

warming_ranking <- all_warm_imp_summary %>%
  dplyr::group_by(TPL,Taxa) %>%
  dplyr::arrange(desc(Mean_RelativeImportance_AUC),.by_group=TRUE) %>%
  dplyr::mutate(Mean_Importance_Rank=row_number()) %>%
  dplyr::ungroup()

write.csv(warming_ranking,
          file.path(all_out,
                    "12_ALL_TAXA_warming_contrast_process_ranking.csv"),
          row.names=FALSE)


## ============================================================
## Cross-taxa figures
## ============================================================

process_labels <- c(
  immigration="Immigration",
  growth="Growth",
  interaction="Interaction",
  noise="Stochastic noise"
)

taxa_labels <- c(
  Bacteria="Bacteria",
  Archaea="Archaea",
  Fungi_18S="Fungi",
  Protists="Protists",
  Nematoda="Nematoda"
)

## Fig 1: condition-specific importance
p1 <- ggplot(
  all_imp_rank,
  aes(x=Process,y=RelativeImportance_AUC)
) +
  geom_boxplot(width=0.58,outlier.shape=NA) +
  geom_jitter(width=0.08,height=0,size=1.2,alpha=0.6) +
  facet_grid(TPL + Condition ~ Taxa,
             labeller=labeller(Taxa=taxa_labels)) +
  scale_x_discrete(limits=process_order,labels=process_labels) +
  labs(x=NULL,y="Relative process-intensity sensitivity") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggsave(
  file.path(all_out,
            "Fig1_ALL_TAXA_condition_specific_relative_process_sensitivity.pdf"),
  p1,width=15,height=8
)

## Fig 2: warming-associated beta-shift importance
p2 <- ggplot(
  all_warm_imp_rank,
  aes(x=Process,y=RelativeImportance_AUC)
) +
  geom_boxplot(width=0.58,outlier.shape=NA) +
  geom_jitter(width=0.08,height=0,size=1.2,alpha=0.6) +
  facet_grid(TPL ~ Taxa,
             labeller=labeller(Taxa=taxa_labels)) +
  scale_x_discrete(limits=process_order,labels=process_labels) +
  labs(
    x=NULL,
    y="Relative sensitivity of warming-associated beta shift"
  ) +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

ggsave(
  file.path(all_out,
            "Fig2_ALL_TAXA_warming_contrast_relative_process_sensitivity.pdf"),
  p2,width=15,height=5.8
)

## Fig 3: warming-associated beta-difference intensity curves
p3 <- ggplot(
  all_warm_rank_curve,
  aes(
    x=Intensity,
    y=DeltaBeta_mean_seed,
    group=interaction(Rank,Process),
    color=Process
  )
) +
  geom_hline(yintercept=0,linetype=2,linewidth=0.35) +
  geom_line(alpha=0.18,linewidth=0.4) +
  geom_line(
    data=all_warm_curve_summary,
    aes(
      x=Intensity,
      y=DeltaBeta_mean_Top10,
      group=Process,
      color=Process
    ),
    inherit.aes=FALSE,
    linewidth=1.05
  ) +
  facet_grid(
    TPL ~ Taxa,
    scales="free_y",
    labeller=labeller(Taxa=taxa_labels)
  ) +
  scale_color_discrete(labels=process_labels) +
  labs(
    x="Process intensity (z)",
    y=expression(Delta*beta==beta[W]-beta[U]),
    color="Process"
  ) +
  theme_classic()

ggsave(
  file.path(all_out,
            "Fig3_ALL_TAXA_warming_contrast_process_intensity_curves.pdf"),
  p3,width=15,height=6.5
)


cat("\n============================================================\n")
cat("ALL taxa tTTPL + pTTPL evaluation finished.\n")
cat("Cross-taxa results saved to:\n",all_out,"\n")
cat("\nMain output tables:\n")
cat("05_ALL_TAXA_process_importance_Top10_summary_for_plot.csv\n")
cat("10_ALL_TAXA_warming_contrast_process_importance_Top10_summary_for_plot.csv\n")
cat("11_ALL_TAXA_condition_specific_process_ranking.csv\n")
cat("12_ALL_TAXA_warming_contrast_process_ranking.csv\n")
cat("============================================================\n")
