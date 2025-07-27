rm(list =ls())
require(ieggr)
##2010to2020
info = lazyopen("./data/summary of environmental data2010-2020.csv")
Plot_list = unique(info$plot)
##Archaea
##Null data
if(!file.exists("Taylor/Archaea/Null/")){
  dir.create("Taylor/Archaea/Null/")
}
beta_summary_all = data.frame(Null = integer(), Plot=character(), Intercept=double(),	Beta=double(), R2=double(), P=double())
nrep=100
for (i in 1:nrep) {
  for(p in 1:length(Plot_list)){
    data_name = paste("./data/Archaea/Null/MPDik_null_", i,"_", Plot_list[p], "_otu_All_abun.csv",sep="")
    subcom = read.csv(data_name, header=T, row.names =1)
    all_summary <- data.frame(Mean = double(), Var = double())
    for(j in 1:nrow(subcom)){
      abun <- as.numeric(subcom[j,])
      mean_abun <- mean(abun)
      var_abun <- var(abun)
      abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
      rownames(abun_summary) <- rownames(subcom)[j]
      all_summary <- rbind(all_summary, abun_summary)
    }
    ##plot real point and fit line
    all_summary <- all_summary[all_summary$Mean>0,]
    all_summary["LogMean"] <- log10(all_summary$Mean)
    all_summary["LogVar"] <- log10(all_summary$Variance)
    model.lm <- lm(LogVar ~ LogMean, data = all_summary)
    l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
              beta = format(coef(model.lm)[[2]]),
              r2 = format(summary(model.lm)$r.squared, digits = 3),
              p = format(summary(model.lm)$coefficients[2,4], digits = 3))
    beta_summary = data.frame(Null = i, Plot=Plot_list[p], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
    beta_summary_all = rbind(beta_summary_all, beta_summary)
  }
}
write.csv(beta_summary_all, "./Taylor/Archaea/Null/MPDik_null_all_otu_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Bacteria
##Null data
if(!file.exists("Taylor/Bacteria/Null/")){
  dir.create("Taylor/Bacteria/Null/")
}
beta_summary_all = data.frame(Null = integer(), Plot=character(), Intercept=double(),	Beta=double(), R2=double(), P=double())
nrep=100
for (i in 1:nrep) {
  for(p in 1:length(Plot_list)){
    data_name = paste("./data/Bacteria/Null/MPDik_null_", i,"_", Plot_list[p], "_otu_All_abun.csv",sep="")
    subcom = read.csv(data_name, header=T, row.names =1)
    all_summary <- data.frame(Mean = double(), Var = double())
    for(j in 1:nrow(subcom)){
      abun <- as.numeric(subcom[j,])
      mean_abun <- mean(abun)
      var_abun <- var(abun)
      abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
      rownames(abun_summary) <- rownames(subcom)[j]
      all_summary <- rbind(all_summary, abun_summary)
    }
    ##plot real point and fit line
    all_summary <- all_summary[all_summary$Mean>0,]
    all_summary["LogMean"] <- log10(all_summary$Mean)
    all_summary["LogVar"] <- log10(all_summary$Variance)
    model.lm <- lm(LogVar ~ LogMean, data = all_summary)
    l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
              beta = format(coef(model.lm)[[2]]),
              r2 = format(summary(model.lm)$r.squared, digits = 3),
              p = format(summary(model.lm)$coefficients[2,4], digits = 3))
    beta_summary = data.frame(Null = i, Plot=Plot_list[p], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
    beta_summary_all = rbind(beta_summary_all, beta_summary)
  }
}
write.csv(beta_summary_all, "./Taylor/Bacteria/Null/MPDik_null_all_otu_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Fungi
##Null data
if(!file.exists("Taylor/Fungi/Null/")){
  dir.create("Taylor/Fungi/Null/")
}
beta_summary_all = data.frame(Null = integer(), Plot=character(), Intercept=double(),	Beta=double(), R2=double(), P=double())
nrep=100
for (i in 1:nrep) {
  for(p in 1:length(Plot_list)){
    data_name = paste("./data/Fungi/Null/MPDik_null_", i,"_", Plot_list[p], "_otu_All_abun.csv",sep="")
    subcom = read.csv(data_name, header=T, row.names =1)
    all_summary <- data.frame(Mean = double(), Var = double())
    for(j in 1:nrow(subcom)){
      abun <- as.numeric(subcom[j,])
      mean_abun <- mean(abun)
      var_abun <- var(abun)
      abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
      rownames(abun_summary) <- rownames(subcom)[j]
      all_summary <- rbind(all_summary, abun_summary)
    }
    ##plot real point and fit line
    all_summary <- all_summary[all_summary$Mean>0,]
    all_summary["LogMean"] <- log10(all_summary$Mean)
    all_summary["LogVar"] <- log10(all_summary$Variance)
    model.lm <- lm(LogVar ~ LogMean, data = all_summary)
    l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
              beta = format(coef(model.lm)[[2]]),
              r2 = format(summary(model.lm)$r.squared, digits = 3),
              p = format(summary(model.lm)$coefficients[2,4], digits = 3))
    beta_summary = data.frame(Null = i, Plot=Plot_list[p], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
    beta_summary_all = rbind(beta_summary_all, beta_summary)
  }
}
write.csv(beta_summary_all, "./Taylor/Fungi/Null/MPDik_null_all_otu_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Protists
##Null data
if(!file.exists("Taylor/Protists/Null/")){
  dir.create("Taylor/Protists/Null/")
}
beta_summary_all = data.frame(Null = integer(), Plot=character(), Intercept=double(),	Beta=double(), R2=double(), P=double())
nrep=100
for (i in 1:nrep) {
  for(p in 1:length(Plot_list)){
    data_name = paste("./data/Protists/Null/MPDik_null_", i,"_", Plot_list[p], "_otu_All_abun.csv",sep="")
    subcom = read.csv(data_name, header=T, row.names =1)
    all_summary <- data.frame(Mean = double(), Var = double())
    for(j in 1:nrow(subcom)){
      abun <- as.numeric(subcom[j,])
      mean_abun <- mean(abun)
      var_abun <- var(abun)
      abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
      rownames(abun_summary) <- rownames(subcom)[j]
      all_summary <- rbind(all_summary, abun_summary)
    }
    ##plot real point and fit line
    all_summary <- all_summary[all_summary$Mean>0,]
    all_summary["LogMean"] <- log10(all_summary$Mean)
    all_summary["LogVar"] <- log10(all_summary$Variance)
    model.lm <- lm(LogVar ~ LogMean, data = all_summary)
    l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
              beta = format(coef(model.lm)[[2]]),
              r2 = format(summary(model.lm)$r.squared, digits = 3),
              p = format(summary(model.lm)$coefficients[2,4], digits = 3))
    beta_summary = data.frame(Null = i, Plot=Plot_list[p], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
    beta_summary_all = rbind(beta_summary_all, beta_summary)
  }
}
write.csv(beta_summary_all, "./Taylor/Protists/Null/MPDik_null_all_otu_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Nematoda
##Null data
if(!file.exists("Taylor/Nematoda/Null/")){
  dir.create("Taylor/Nematoda/Null/")
}
beta_summary_all = data.frame(Null = integer(), Plot=character(), Intercept=double(),	Beta=double(), R2=double(), P=double())
nrep=100
for (i in 1:nrep) {
  for(p in 1:length(Plot_list)){
    data_name = paste("./data/Nematoda/Null/MPDik_null_", i,"_", Plot_list[p], "_otu_All_abun.csv",sep="")
    subcom = read.csv(data_name, header=T, row.names =1)
    all_summary <- data.frame(Mean = double(), Var = double())
    for(j in 1:nrow(subcom)){
      abun <- as.numeric(subcom[j,])
      mean_abun <- mean(abun)
      var_abun <- var(abun)
      abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
      rownames(abun_summary) <- rownames(subcom)[j]
      all_summary <- rbind(all_summary, abun_summary)
    }
    ##plot real point and fit line
    all_summary <- all_summary[all_summary$Mean>0,]
    all_summary["LogMean"] <- log10(all_summary$Mean)
    all_summary["LogVar"] <- log10(all_summary$Variance)
    model.lm <- lm(LogVar ~ LogMean, data = all_summary)
    l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
              beta = format(coef(model.lm)[[2]]),
              r2 = format(summary(model.lm)$r.squared, digits = 3),
              p = format(summary(model.lm)$coefficients[2,4], digits = 3))
    beta_summary = data.frame(Null = i, Plot=Plot_list[p], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
    beta_summary_all = rbind(beta_summary_all, beta_summary)
  }
}
write.csv(beta_summary_all, "./Taylor/Nematoda/Null/MPDik_null_all_otu_Taylor_beta.csv", quote = FALSE, row.names = FALSE)
