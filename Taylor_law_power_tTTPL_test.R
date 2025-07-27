rm(list=ls())
##Test the Taylor law(all ASVs/OTUs)
require(ggplot2)
require(ieggr)

if(!file.exists("Taylor/")){
  dir.create("Taylor")
}
##Time:Type IV PLE for mixed-species population temporal aggregation (stability)
if(!file.exists("Taylor/")){
  dir.create("Taylor/")
}
if(!file.exists("Taylor/Archaea/")){
  dir.create("Taylor/Archaea/")
}
if(!file.exists("Taylor/Bacteria/")){
  dir.create("Taylor/Bacteria/")
}
if(!file.exists("Taylor/Protists/")){
  dir.create("Taylor/Protists/")
}
if(!file.exists("Taylor/Nematoda/")){
  dir.create("Taylor/Nematoda/")
}
if(!file.exists("Taylor/Fungi/")){
  dir.create("Taylor/Fungi/")
}

info = lazyopen("./data/summary of environmental data2010-2020.csv")
Plot_list = unique(info$plot)
##Archaea
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Archaea/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Archaea/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Archaea/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Archaea/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)

##Bacteria
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Bacteria/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Bacteria/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Bacteria/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Bacteria/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)

##Fungi
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Fungi/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Fungi/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Fungi/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Fungi/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)
##Fungi_18S
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Fungi_18S/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Fungi_18S/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Fungi_18S/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Fungi_18S/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)

##Protists
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Protists/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Protists/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Protists/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Protists/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)

##Nematoda
beta_summary_all = data.frame(Plot=character(),	Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Nematoda/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
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
  all_summary["LogMean"] <- log10(all_summary$Mean)
  all_summary["LogVar"] <- log10(all_summary$Variance)
  summary_name <- paste("./Taylor/Nematoda/",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
  write.csv(all_summary, summary_name, quote = FALSE)
  plot_name <- paste("./Taylor/Nematoda/",Plot_list[i],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
  par(mar=c(3,5,2,2))
  p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
    geom_point(size=1)  +
    geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
  model.lm <- lm(LogVar ~ LogMean, data = all_summary)
  l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
            beta = format(coef(model.lm)[[2]]),
            r2 = format(summary(model.lm)$r.squared, digits = 3),
            p = format(summary(model.lm)$coefficients[2,4], digits = 3))
  eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
  p + theme_bw() +
    geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
    labs(x="Log mean abundance", y = "Log variance")
  dev.off()
  ggsave(plot_name, width = 4, height = 4)
  dev.new()
  beta_summary = data.frame(Plot=Plot_list[i], Intercept=as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
  beta_summary_all = rbind(beta_summary_all, beta_summary)
}
write.csv(beta_summary_all, "./Taylor/Nematoda/UPARSE_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)