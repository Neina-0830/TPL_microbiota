# Standardized mixed-effects associations. Ecosystem function is the response.
source(file.path(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),"code/analysis_helpers.R"))
setwd(project_root());require_packages(c("lmerTest","pbkrtest"))
x <- read.csv("results/10_environment_SEM_inputs/all_environment_TPL_function_48plots.csv",row.names=1,check.names=FALSE)
x$Plot <- rownames(x);x <- design_factors(x)
map <- data.frame(Group=rep(c("Bacteria","Archaea","Fungi","Protists","Nematoda"),2),
 TPL=rep(c("tTTPL","pTTPL"),each=5),Beta=c("Bac_beta","Arch_beta","Fun_18S_beta","Pro_beta","Nem_beta",
 "Bac_pTTPL_beta","Arch_pTTPL_beta","Fun_18S_pTTPL_beta","Pro_pTTPL_beta","Nem_pTTPL_beta"))
functions <- unlist(lapply(c("GPP","NEE","ER","Ra","Rh","Rs"),function(z) paste(z,c("mean","cv"),sep="_")))
envs <- grep("^(temperature|moisture|pH|TC|NO3N|NH4N|plant_biomass|plant_richness)_(mean|cv)$",names(x),value=TRUE)
res <- list()
for(i in seq_len(nrow(map))) for(kind in c("environment","function")) {
  vars <- if(kind=="environment") envs else functions
  for(v in vars) {
    response <- if(kind=="environment") map$Beta[i] else v
    predictor <- if(kind=="environment") v else map$Beta[i]
    d <- x[complete.cases(x[,c(response,predictor)]),,drop=FALSE]
    if(nrow(d)<20 || sd(d[[response]])==0 || sd(d[[predictor]])==0) stop("Insufficient variation: ",response," / ",predictor)
    d$y <- as.numeric(scale(d[[response]]));d$z <- as.numeric(scale(d[[predictor]]))
    fit <- lmerTest::lmer(y~z+(1|Block/WholePlotID),d)
    tab <- coef(summary(fit,ddf="Kenward-Roger"))["z",]
    aa <- anova(fit,type=2,ddf="Kenward-Roger")
    res[[paste(i,kind,v)]] <- data.frame(Group=map$Group[i],TPL=map$TPL[i],Kind=kind,
      Variable=v,Response=response,Predictor=predictor,N=nrow(d),b=tab["Estimate"],
      SE=tab["Std. Error"],df=tab["df"],F=aa$`F value`[1],P=aa$`Pr(>F)`[1],
      CI_low=tab["Estimate"]-qt(.975,tab["df"])*tab["Std. Error"],
      CI_high=tab["Estimate"]+qt(.975,tab["df"])*tab["Std. Error"])
  }
}
res <- do.call(rbind,res)
# Raw P values are preserved; Holm values are supplied by group, type and family.
res$P_Holm <- ave(res$P,interaction(res$Group,res$TPL,res$Kind),FUN=function(p) p.adjust(p,"holm"))
write_table(res,"results/13_ecosystem_function/mixed_effects_associations.csv")
