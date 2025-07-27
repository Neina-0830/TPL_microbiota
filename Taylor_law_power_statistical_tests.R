rm(list=ls())
require(ggplot2)
require(ieggr)
library(dplyr)
require(reshape2)
data_summary <-  function(data=NULL, varname, groupnames=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  #group_by + summarise
  datac <- ddply(data, groupnames, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 varname
  )
  datac <- plyr::rename(datac, c("mean" = varname))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

if(!file.exists("Taylor/compare/")){
  dir.create("Taylor/compare/")
}
if(!file.exists("Taylor/Archaea/compare/")){
  dir.create("Taylor/Archaea/compare/")
}
if(!file.exists("Taylor/Bacteria/compare/")){
  dir.create("Taylor/Bacteria/compare/")
}
if(!file.exists("Taylor/Protists/compare/")){
  dir.create("Taylor/Protists/compare/")
}
if(!file.exists("Taylor/Nematoda/compare/")){
  dir.create("Taylor/Nematoda/compare/")
}
if(!file.exists("Taylor/Fungi/compare/")){
  dir.create("Taylor/Fungi/compare/")
}

Plot_info = lazyopen("./data/treatment_list10_20.txt")
Plot_info$Plot = rownames(Plot_info)
##tTTPL
##Archaea
Taylor_para = read.csv("./Taylor/Archaea/UPARSE_all_otus_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Archaea/UPARSE_all_otus_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Archaea/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Bacteria
Taylor_para = read.csv("./Taylor/Bacteria/UPARSE_all_otus_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Bacteria/UPARSE_all_otus_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Bacteria/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Fungi
Taylor_para = read.csv("./Taylor/Fungi/UPARSE_all_otus_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Fungi/UPARSE_all_otus_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Fungi/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Protists
Taylor_para = read.csv("./Taylor/Protists/UPARSE_all_otus_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Protists/UPARSE_all_otus_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Protists/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Nematoda
Taylor_para = read.csv("./Taylor/Nematoda/UPARSE_all_otus_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Nematoda/UPARSE_all_otus_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Nematoda/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##pTTPL
##Archaea
Taylor_para = read.csv("./Taylor/Archaea/MPD_weighted_otu_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Archaea/MPD_weighted_otu_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Archaea/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Bacteria
Taylor_para = read.csv("./Taylor/Bacteria/MPD_weighted_otu_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Bacteria/MPD_weighted_otu_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Bacteria/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Fungi
Taylor_para = read.csv("./Taylor/Fungi/MPD_weighted_otu_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Fungi/MPD_weighted_otu_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Fungi/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Protists
Taylor_para = read.csv("./Taylor/Protists/MPD_weighted_otu_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Protists/MPD_weighted_otu_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Protists/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)


##Nematoda
Taylor_para = read.csv("./Taylor/Nematoda/MPD_weighted_otu_Taylor_beta.csv")
Taylor_para_info = merge(Taylor_para, Plot_info, by="Plot")
write.csv(Taylor_para_info, "./Taylor/Nematoda/MPD_weighted_otu_Taylor_beta_info.csv", quote = F)
##statistic analysis
Treatment_list = unique(Taylor_para_info$Treatment)
compre_result_all = as.data.frame(matrix(nrow=0, ncol=27))
for(i in 1:length(Treatment_list)){
  Taylor_para_info_control = Taylor_para_info[Taylor_para_info$Treatment=="Ctrl",]
  Taylor_para_info_treatment = Taylor_para_info[Taylor_para_info$Treatment==Treatment_list[i],]
  compre_result = as.data.frame(compv(Taylor_para_info_treatment$Beta, Taylor_para_info_control$Beta, paired = FALSE, alternative = c("two.sided", "less", "greater"),
                                      method = c("t", "w", "d", "r"), nworker = 1, rr.fill.zero = c("yes", "no")))
  compre_result["Treatment.x"] = Treatment_list[i]
  compre_result["Treatment.y"] = "Ctrl"
  compre_result_all = rbind(compre_result_all,compre_result)
}
write.csv(compre_result_all, "./Taylor/Nematoda/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", quote = F, row.names = F)