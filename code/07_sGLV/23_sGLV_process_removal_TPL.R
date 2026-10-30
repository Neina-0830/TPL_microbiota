# Fit tTTPL and pTTPL for the full model and four process removals.
source(file.path(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),"code/analysis_helpers.R"))
setwd(project_root())
root <- normalizePath(Sys.getenv("TPL_SGLV_WORK_DIR",unset="results/17_sGLV"),mustWork=FALSE)
cfg <- read.csv("code/07_sGLV/sGLV_group_config.csv",check.names=FALSE)
groups <- Sys.getenv("TPL_SGLV_GROUPS",unset="all")
groups <- if(groups=="all") cfg$Group else strsplit(groups,",",fixed=TRUE)[[1]]
scenarios <- c("full","immigration_removed","growth_removed","interaction_removed","noise_removed")
tpls <- strsplit(Sys.getenv("TPL_SGLV_TPL",unset="tTTPL,pTTPL"),",",fixed=TRUE)[[1]]
all <- list()
for(group in groups) for(tpl in tpls) {
  if(tpl=="pTTPL") {
    require_packages(c("iCAMP","bigmemory"))
    pdfile <- file.path("results/pairwise_distance",if(group=="Fungi_18S") "Fungi" else group,
                        "All_OTUs_2010to2020_pdist_out.rda")
    e <- new.env();load(pdfile,envir=e)
    pd <- e$pdist.out
    if(is.null(pd)) stop("Missing pdist.out in ",pdfile)
    pd$pd.wd <- normalizePath(dirname(pdfile),mustWork=TRUE)
    desc <- file.path(pd$pd.wd,pd$pd.file)
    big <- tryCatch(bigmemory::attach.big.matrix(dget(desc)),
                    error=function(e) bigmemory::attach.big.matrix(desc))
  }
  ans <- list()
  for(cond in c("warming","unwarming")) for(scenario in scenarios) {
    folder <- file.path(root,"Resample",group,tpl,"process_removal_top10",scenario,cond)
    files <- list.files(folder,pattern="\\.csv$",full.names=TRUE)
    if(length(files)!=30L) stop("Expected 10 parameter sets x 3 seeds in ",folder,"; found ",length(files))
    for(file in files) {
      name <- basename(file)
      hit <- regmatches(name,regexec("^Time_series_from_sGLV_rank([0-9]+)_.*_seed([0-9]+)_.*\\.csv$",name))[[1]]
      if(length(hit)!=3) stop("Cannot parse ",name)
      x <- as.matrix(read.csv(file,row.names=1,check.names=FALSE));storage.mode(x) <- "double"
      if(!all(paste0("Y",2010:2020) %in% colnames(x))) stop("Missing simulation years: ",file)
      x <- x[,paste0("Y",2010:2020),drop=FALSE]
      if(any(!is.finite(x)) || any(x<0) || any(colSums(x)<=0)) stop("Invalid simulation: ",file)
      x <- sweep(x,2,colSums(x),"/")
      if(tpl=="pTTPL") {
        ids <- match(rownames(x),pd$tip.label)
        if(anyNA(ids)) stop("Simulated taxa are absent from the phylogeny.")
        # Same abundance-weighted MPD normalization as the empirical module.
        contributions <- matrix(0,nrow(x),ncol(x),dimnames=dimnames(x))
        chunk <- 500L
        for(start in seq(1L,nrow(x),by=chunk)) {
          ii <- start:min(start+chunk-1L,nrow(x))
          distances <- big[ids[ii],ids,drop=FALSE]
          contributions[ii,] <- x[ii,,drop=FALSE]*(distances%*%x)
        }
        total <- colSums(contributions)
        if(any(total<=0)) stop("Zero abundance-weighted phylogenetic distance: ",file)
        x <- sweep(contributions,2,total,"/")
      }
      fit <- tpl_fit(rowMeans(x),apply(x,1,var))
      ans[[paste(cond,scenario,hit[2],hit[3])]] <- data.frame(Group=group,TPL=tpl,
        Condition=cond,Scenario=scenario,Rank=as.integer(hit[2]),Seed=as.integer(hit[3]),
        Beta=unname(fit["Beta"]),SE_Beta=unname(fit["SE_Beta"]),R2=unname(fit["R2"]))
    }
  }
  a <- do.call(rbind,ans)
  if(anyDuplicated(a[c("Group","TPL","Condition","Scenario","Rank","Seed")])) stop("Duplicate simulation keys.")
  write_table(a,file.path(root,"Taylor",group,tpl,"process_removal_top10","TPL_estimates.csv"))
  all[[paste(group,tpl)]] <- a
}
write_table(do.call(rbind,all),file.path(root,"Taylor/process_removal_top10/TPL_estimates.csv"))
