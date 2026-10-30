# Temporal variance-mean scaling, detrending and taxon-removal robustness.
source(file.path(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),"code/analysis_helpers.R"))
setwd(project_root())
meta <- read.csv("data/prepared/environment_metadata.csv",check.names=FALSE)
out <- "results/01_tTTPL";dir.create(out,recursive=TRUE,showWarnings=FALSE)
groups <- c("Bacteria","Archaea","Fungi","Protists","Nematoda")
reps <- as.integer(Sys.getenv("TPL_SUBCOMPOSITION_REPS",unset="200"))
set.seed(20260813)
all <- sensitivity <- concentration <- references <- list()
for(g in groups) {
  plots <- unique(meta$plot)
  mats <- lapply(plots,function(p) {
    md <- meta[meta$plot==p & meta$year %in% 2010:2020,,drop=FALSE];md <- md[order(md$year),]
    if(nrow(md)!=11 || anyDuplicated(md$year)) stop("Expected 11 annual samples: ",g," ",p)
    f <- file.path("data/prepared/plot_abundance",g,paste0(p,"_UPARSE_otu_All_abun.csv"))
    x <- as.matrix(read.csv(f,row.names=1,check.names=FALSE));storage.mode(x) <- "double"
    x <- x[,md$Sample_name,drop=FALSE]
    if(any(!is.finite(x)) || any(x<0)) stop("Invalid abundance values: ",f)
    list(mat=x,md=md)
  });names(mats) <- plots
  pooled <- unlist(lapply(mats,function(z) {m <- rowMeans(z$mat);log10(m[m>0])}))
  ref <- median(pooled);references[[g]] <- data.frame(Group=g,Log10MeanReference=ref,MeanAbundanceReference=10^ref)
  # Group-wide taxon ranking by mean absolute temporal trend across subplots.
  taxa <- unique(unlist(lapply(mats,function(z) rownames(z$mat))))
  trends <- matrix(0,length(taxa),length(plots),dimnames=list(taxa,plots))
  for(p in plots) {
    z <- mats[[p]];t <- z$md$year-mean(z$md$year)
    slopes <- drop(z$mat%*%t/sum(t^2));trends[rownames(z$mat),p] <- abs(slopes)
    score <- sort(abs(slopes),decreasing=TRUE)
    concentration[[paste(g,p)]] <- data.frame(Group=g,Plot=p,TopK=c(1,5,10),
      Trend_fraction=vapply(c(1,5,10),function(k) sum(head(score,k))/sum(score),numeric(1)))
  }
  ranking <- names(sort(rowMeans(trends),decreasing=TRUE))
  subsets <- list()
  for(k in c(1,5)) subsets[[paste0("targeted_top",k)]] <- head(ranking,k)
  for(i in seq_len(min(10,length(ranking)))) subsets[[paste0("leave_one_out_",i)]] <- ranking[i]
  for(frac in c(.01,.05,.10)) for(r in seq_len(reps)) {
    subsets[[paste0("random_",frac,"_",r)]] <- sample(taxa,max(1L,ceiling(length(taxa)*frac)))
  }
  for(p in plots) {
    z <- mats[[p]];x <- z$mat;md <- z$md;t <- md$year-mean(md$year)
    mu <- rowMeans(x);va <- apply(x,1,var)
    residuals <- x-mu-outer(drop(x%*%t/sum(t^2)),t)
    raw <- tpl_fit(mu,va,ref);det <- tpl_fit(mu,apply(residuals,1,var),ref)
    d <- data.frame(Group=g,Plot=p,Log10MeanReference=ref)
    for(nm in names(raw)) d[[paste0(nm,"_raw")]] <- raw[nm]
    for(nm in names(det)) d[[paste0(nm,"_detrended")]] <- det[nm]
    for(nm in c("block","Number","NS","warm","precip","clip","Treatment","WholePlotID","Precipitation")) if(nm %in% names(md)) d[[nm]] <- md[[nm]][1]
    all[[paste(g,p)]] <- d
    write_table(data.frame(Taxon=rownames(x),Mean=mu,Variance=va,Detrended_variance=apply(residuals,1,var)),
                file.path(out,"taxon_level_tables",paste0(g,"_",p,".csv")))
    for(nm in names(subsets)) {
      keep <- !rownames(x) %in% subsets[[nm]];a <- x[keep,,drop=FALSE]
      if(any(colSums(a)<=0)) stop("Taxon removal left an empty annual community.")
      a <- sweep(a,2,colSums(a),"/")
      beta <- tpl_fit(rowMeans(a),apply(a,1,var),ref)["Beta"]
      sensitivity[[paste(g,p,nm)]] <- data.frame(Group=g,Plot=p,Scenario=nm,Beta=unname(beta))
    }
  }
}
x <- do.call(rbind,all)
write_table(x,file.path(out,"UPARSE_all_groups_plot_level_revision_metrics.csv"))
saveRDS(x,file.path(out,"UPARSE_all_groups_plot_level_revision_metrics.rds"))
write_table(do.call(rbind,references),file.path(out,"group_abundance_reference.csv"))
write_table(do.call(rbind,concentration),file.path(out,"trend_concentration.csv"))
sens <- do.call(rbind,sensitivity)
write_table(sens,file.path(out,"taxon_removal_plot_betas.csv"))
require_packages(c("lmerTest","emmeans","pbkrtest"))
design <- read.csv("data/prepared/treatment_metadata_48plots.csv",check.names=FALSE)
sens <- merge(sens,design,by.x="Plot",by.y="plot",all.x=TRUE);sens <- design_factors(sens)
ans <- lapply(split(sens,interaction(sens$Group,sens$Scenario,drop=TRUE)),function(d) {
  fit <- lmerTest::lmer(Beta~WarmF*PrecipF*ClipF+(1|Block/WholePlotID),d)
  a <- overall_contrasts(fit);cbind(Group=d$Group[1],Scenario=d$Scenario[1],a[a$contrast=="Warming",])
})
write_table(do.call(rbind,ans),file.path(out,"taxon_removal_warming_effects.csv"))
