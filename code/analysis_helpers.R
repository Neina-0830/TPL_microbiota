# Shared fitting and split-plot inference used by the empirical modules.
project_root <- function() normalizePath(Sys.getenv("TPL_PROJECT_ROOT", unset=getwd()), mustWork=TRUE)
write_table <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  old <- options(digits=17); on.exit(options(old))
  write.csv(x, path, row.names=FALSE, quote=TRUE, na="")
}
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
  if(length(missing)) stop("Install required packages: ", paste(missing, collapse=", "))
}
design_factors <- function(d) {
  d$Block <- factor(d$block)
  if(!"WholePlotID" %in% names(d)) d$WholePlotID <- paste(d$block, sub("^[NS]", "", d$Plot), sep="_")
  d$WholePlotID <- factor(d$WholePlotID)
  d$WarmF <- factor(d$warm, levels=c(0,1), labels=c("Unwarmed","Warmed"))
  if("Precipitation" %in% names(d)) {
    p <- tolower(as.character(d$Precipitation))
    p[p %in% c("normal","ambient")] <- "Ambient"
    p[p %in% c("half","reduced")] <- "Half"
    p[p %in% c("double","increased")] <- "Double"
  } else {
    p <- as.character(d$PrecipF)
  }
  d$PrecipF <- factor(p, levels=c("Ambient","Half","Double"))
  d$ClipF <- factor(d$clip, levels=c(0,1), labels=c("Unclipped","Clipped"))
  if(anyNA(d[c("Block","WholePlotID","WarmF","PrecipF","ClipF")])) stop("Incomplete experimental design.")
  d
}
tpl_fit <- function(mu, va, reference=0) {
  keep <- is.finite(mu) & is.finite(va) & mu>0 & va>0
  if(sum(keep)<3L || sd(log10(mu[keep]))==0) stop("Insufficient taxa for a finite TPL fit.")
  d <- data.frame(x=log10(mu[keep])-reference, y=log10(va[keep]))
  fit <- lm(y~x,d); s <- summary(fit); v <- vcov(fit)
  c(Intercept=unname(coef(fit)[1]-reference*coef(fit)[2]),
    Beta=unname(coef(fit)[2]), SE_Beta=s$coefficients[2,2], R2=s$r.squared,
    P_Beta=s$coefficients[2,4], LogVar_at_ref=unname(coef(fit)[1]),
    SE_LogVar_at_ref=s$coefficients[1,2], Cov_LogVarRef_Beta=v[1,2], n_taxa=sum(keep))
}
overall_contrasts <- function(fit) {
  a <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit,~WarmF,weights="equal"),list(Warming=c(-1,1))))
  b <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit,~PrecipF,weights="equal"),
       list(Reduced_precipitation=c(-1,1,0),Increased_precipitation=c(-1,0,1))))
  c <- as.data.frame(emmeans::contrast(emmeans::emmeans(fit,~ClipF,weights="equal"),list(Clipping=c(-1,1))))
  x <- rbind(a,b,c); x$P_Holm <- p.adjust(x$p.value,"holm"); x
}
treatment_contrasts <- function(fit) {
  emm <- emmeans::emmeans(fit,~WarmF*PrecipF*ClipF)
  grid <- as.data.frame(emm)
  control <- which(grid$WarmF=="Unwarmed" & grid$PrecipF=="Ambient" & grid$ClipF=="Unclipped")
  x <- as.data.frame(emmeans::contrast(emm, method="trt.vs.ctrl",ref=control,adjust="none"))
  x$P_Holm <- p.adjust(x$p.value,"holm"); x
}
run_treatment_statistics <- function(tpl) {
  require_packages(c("lme4","lmerTest","emmeans","pbkrtest","metafor"))
  root <- project_root(); setwd(root)
  file <- if(tpl=="tTTPL") "results/01_tTTPL/UPARSE_all_groups_plot_level_revision_metrics.csv" else "results/04_pTTPL/pTTPL_all_groups_plot_level_full_precision.csv"
  x <- read.csv(file,check.names=FALSE)
  if(!"Plot" %in% names(x)) x$Plot <- x$plot
  meta <- read.csv("data/prepared/treatment_metadata_48plots.csv",check.names=FALSE)
  # Treatment metadata are authoritative and joined by subplot identifier.
  cols <- setdiff(names(meta),"plot")
  x <- x[,setdiff(names(x),cols),drop=FALSE]
  x <- merge(x,meta,by.x="Plot",by.y="plot",all.x=TRUE,sort=FALSE)
  x <- design_factors(x)
  out <- paste0("results/",if(tpl=="tTTPL") "02_tTTPL" else "05_pTTPL","_factorial_statistics")
  overall <- specific <- tests <- uncertainty <- effects <- list()
  for(g in unique(x$Group)) {
    d <- x[x$Group==g,,drop=FALSE]
    metrics <- intersect(c("Beta_raw","Beta_detrended","LogVar_at_ref_raw","LogVar_at_ref_detrended"),names(d))
    for(metric in metrics) {
      dd <- d[is.finite(d[[metric]]),,drop=FALSE]; dd$y <- dd[[metric]]
      fit <- lmerTest::lmer(y~WarmF*PrecipF*ClipF+(1|Block/WholePlotID),dd,REML=TRUE)
      tag <- paste(g,metric,sep="__")
      overall[[tag]] <- cbind(Group=g,Metric=metric,overall_contrasts(fit))
      specific[[tag]] <- cbind(Group=g,Metric=metric,treatment_contrasts(fit))
      aa <- as.data.frame(anova(fit,type=2,ddf="Kenward-Roger"));aa$Term <- rownames(aa)
      tests[[tag]] <- cbind(Group=g,Metric=metric,aa)
      dir.create(file.path(out,"models"),recursive=TRUE,showWarnings=FALSE)
      saveRDS(fit,file.path(out,"models",paste0(tag,".rds")))
      # Descriptive pooled-SD effect sizes for the same four marginal factors.
      selectors <- list(Warming=list(dd$warm==1,dd$warm==0),
        Reduced_precipitation=list(dd$PrecipF=="Half",dd$PrecipF=="Ambient"),
        Increased_precipitation=list(dd$PrecipF=="Double",dd$PrecipF=="Ambient"),
        Clipping=list(dd$clip==1,dd$clip==0))
      effects[[tag]] <- do.call(rbind,lapply(names(selectors),function(nm) {
        a <- dd$y[selectors[[nm]][[1]]];b <- dd$y[selectors[[nm]][[2]]]
        sp <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
        data.frame(Group=g,Metric=metric,Effect=nm,Cohen_d=(mean(a)-mean(b))/sp)
      }))
    }
    # Propagate slope-estimation uncertainty for overall treatment contrasts only.
    dd <- d[is.finite(d$Beta_raw)&is.finite(d$SE_Beta_raw)&d$SE_Beta_raw>0,,drop=FALSE]
    dd$BlockWholePlotID <- interaction(dd$Block,dd$WholePlotID,drop=TRUE)
    mm <- metafor::rma.mv(yi=Beta_raw,V=SE_Beta_raw^2,mods=~WarmF*PrecipF*ClipF,
      random=list(~1|Block,~1|BlockWholePlotID,~1|Plot),data=dd,method="REML",test="t",dfs="contain")
    grid <- expand.grid(WarmF=levels(dd$WarmF),PrecipF=levels(dd$PrecipF),ClipF=levels(dd$ClipF))
    for(nm in c("WarmF","PrecipF","ClipF")) grid[[nm]] <- factor(grid[[nm]],levels=levels(dd[[nm]]))
    M <- model.matrix(~WarmF*PrecipF*ClipF,grid)
    rows <- list(Warming=list(grid$WarmF=="Warmed",grid$WarmF=="Unwarmed"),
      Reduced_precipitation=list(grid$PrecipF=="Half",grid$PrecipF=="Ambient"),
      Increased_precipitation=list(grid$PrecipF=="Double",grid$PrecipF=="Ambient"),
      Clipping=list(grid$ClipF=="Clipped",grid$ClipF=="Unclipped"))
    z <- do.call(rbind,lapply(names(rows),function(nm) {
      L <- colMeans(M[rows[[nm]][[1]],,drop=FALSE])-colMeans(M[rows[[nm]][[2]],,drop=FALSE])
      aa <- anova(mm,X=matrix(L,nrow=1)); data.frame(Group=g,Effect=nm,
        estimate=sum(L*coef(mm)),SE=sqrt(drop(L%*%vcov(mm)%*%L)),P=aa$pval[1])
    }));z$P_Holm <- p.adjust(z$P,"holm");uncertainty[[g]] <- z
  }
  for(nm in c("overall","specific","tests","uncertainty","effects")) write_table(do.call(rbind,get(nm)),file.path(out,paste0(nm,".csv")))
  require_packages("ggplot2")
  est <- do.call(rbind,overall)
  ef <- do.call(rbind,effects)
  a <- merge(ef,est,by.x=c("Group","Metric","Effect"),by.y=c("Group","Metric","contrast"))
  a$Significance <- ifelse(a$P_Holm<.001,"***",ifelse(a$P_Holm<.01,"**",ifelse(a$P_Holm<.05,"*","")))
  write_table(a,file.path(out,"overall_effect_sizes_and_inference.csv"))
  p <- ggplot2::ggplot(a[a$Metric=="Beta_raw",],ggplot2::aes(Effect,Group,fill=Cohen_d))+
    ggplot2::geom_tile()+ggplot2::geom_text(ggplot2::aes(label=Significance))+
    ggplot2::scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B")+
    ggplot2::theme_bw()+ggplot2::labs(x=NULL,y=NULL,fill="Cohen's d",title=tpl)
  ggplot2::ggsave(file.path(out,"overall_treatment_effects.pdf"),p,width=7,height=4)
  # Within-block warming contrasts in beta and fitted variance at a common mean.
  if("LogVar_at_ref_raw" %in% names(x)) {
    av <- aggregate(cbind(Beta_raw,LogVar_at_ref_raw)~Group+Block+PrecipF+ClipF+WarmF,x,mean)
    w <- av[av$WarmF=="Warmed",];u <- av[av$WarmF=="Unwarmed",]
    both <- merge(w,u,by=c("Group","Block","PrecipF","ClipF"),suffixes=c("_w","_u"))
    both$Delta_Beta <- both$Beta_raw_w-both$Beta_raw_u
    both$Delta_LogVar <- both$LogVar_at_ref_raw_w-both$LogVar_at_ref_raw_u
    write_table(both,file.path(out,"warming_beta_reference_variance_changes.csv"))
    assoc <- lapply(split(both,both$Group),function(d) {
      fit <- lmerTest::lmer(Delta_LogVar~scale(Delta_Beta)+PrecipF*ClipF+(1|Block),d)
      a <- as.data.frame(coef(summary(fit,ddf="Kenward-Roger")))
      data.frame(Group=d$Group[1],Term=rownames(a),a,check.names=FALSE)
    })
    write_table(do.call(rbind,assoc),file.path(out,"beta_reference_variance_associations.csv"))
  }
}
