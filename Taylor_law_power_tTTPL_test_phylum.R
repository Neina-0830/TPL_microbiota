rm(list=ls())
##Test the Taylor law(group by different taxon)
require(ggplot2)
require(ieggr)

if(!file.exists("Taylor/Archaea/Phylum/")){
  dir.create("Taylor/Archaea/Phylum/")
}
if(!file.exists("Taylor/Bacteria/Phylum/")){
  dir.create("Taylor/Bacteria/Phylum/")
}
if(!file.exists("Taylor/Fungi/Phylum/")){
  dir.create("Taylor/Fungi/Phylum/")
}
if(!file.exists("Taylor/Protists/Phylum/")){
  dir.create("Taylor/Protists/Phylum/")
}
if(!file.exists("Taylor/Nematoda/Trophic/")){
  dir.create("Taylor/Nematoda/Trophic/")
}

info = lazyopen("./data/summary of environmental data2010-2020.csv")
Plot_list = unique(info$plot)
##Archaea data
##UPARSE
##Group by taxonomy(Phylum level)
beta_summary_all = data.frame(Plot=character(),	Phylum=character(), Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Archaea/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
  if(file.exists(data_name)){
    subcom = read.csv(data_name, header=T, row.names =1)
    taxonomy = read.csv(paste("./data/Archaea/", Plot_list[i], "_UPARSE_otu_All_taxonomy_abun_freq.csv",sep=""), header=T, row.names =1)
    Phylum_num = as.data.frame(table(taxonomy$Phylum))
    colnames(Phylum_num) = c("Phylum","Number")
    write.csv(Phylum_num, paste("./Taylor/Archaea/Phylum/", Plot_list[i], "_UPARSE_otu_All_taxonomy_Phylum_number.csv",sep=""), quote = F, row.names = F)
    Phylum_list = as.character(Phylum_num[Phylum_num$Number>=15,]$Phylum)
    for(phy in 1:length(Phylum_list)){
      subcom_select = subcom[rownames(taxonomy[taxonomy$Phylum==Phylum_list[phy],]),]
      all_summary <- data.frame(Mean = double(), Var = double())
      for(n in 1:nrow(subcom_select)){
        abun <- as.numeric(subcom_select[n,])
        mean_abun <- mean(abun)
        var_abun <- var(abun)
        abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
        rownames(abun_summary) <- rownames(subcom_select)[n]
        all_summary <- rbind(all_summary, abun_summary)
      }
      ##plot real point and fit line
      all_summary = all_summary[all_summary$Variance>0,]
      all_summary["LogMean"] <- log10(all_summary$Mean)
      all_summary["LogVar"] <- log10(all_summary$Variance)
      summary_name <- paste("./Taylor/Archaea/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
      write.csv(all_summary, summary_name, quote = FALSE)
      # plot_name <- paste("./Taylor/Archaea/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
      # par(mar=c(3,5,2,2))
      # p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
      #   geom_point(size=1)  +
      #   geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
      model.lm <- lm(LogVar ~ LogMean, data = all_summary)
      l <- list(beta = format(coef(model.lm)[[2]]),
                r2 = format(summary(model.lm)$r.squared, digits = 3),
                p = format(summary(model.lm)$coefficients[2,4], digits = 3))
      # eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
      # p + theme_bw() +
      #   geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
      #   labs(x="Log mean abundance", y = "Log variance")
      # dev.off()
      # ggsave(plot_name, width = 4, height = 4)
      # dev.new()
      beta_summary = data.frame(Plot=Plot_list[i], Phylum=Phylum_list[phy], Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
      beta_summary_all = rbind(beta_summary_all, beta_summary)
    }
  }
}
write.csv(beta_summary_all, "./Taylor/Archaea/Phylum/UPARSE_diff_Phylum_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Bacteria data
##UPARSE
##Group by taxonomy(Phylum level)
beta_summary_all = data.frame(Plot=character(),	Phylum=character(), Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Bacteria/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
  if(file.exists(data_name)){
    subcom = read.csv(data_name, header=T, row.names =1)
    taxonomy = read.csv(paste("./data/Bacteria/", Plot_list[i], "_UPARSE_otu_All_taxonomy_abun_freq.csv",sep=""), header=T, row.names =1)
    Phylum_num = as.data.frame(table(taxonomy$Phylum))
    colnames(Phylum_num) = c("Phylum","Number")
    write.csv(Phylum_num, paste("./Taylor/Bacteria/Phylum/", Plot_list[i], "_UPARSE_otu_All_taxonomy_Phylum_number.csv",sep=""), quote = F, row.names = F)
    Phylum_list = as.character(Phylum_num[Phylum_num$Number>15,]$Phylum)
    for(phy in 1:length(Phylum_list)){
      subcom_select = subcom[rownames(taxonomy[taxonomy$Phylum==Phylum_list[phy],]),]
      all_summary <- data.frame(Mean = double(), Var = double())
      for(n in 1:nrow(subcom_select)){
        abun <- as.numeric(subcom_select[n,])
        mean_abun <- mean(abun)
        var_abun <- var(abun)
        abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
        rownames(abun_summary) <- rownames(subcom_select)[n]
        all_summary <- rbind(all_summary, abun_summary)
      }
      ##plot real point and fit line
      all_summary = all_summary[all_summary$Variance>0,]
      all_summary["LogMean"] <- log10(all_summary$Mean)
      all_summary["LogVar"] <- log10(all_summary$Variance)
      summary_name <- paste("./Taylor/Bacteria/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
      write.csv(all_summary, summary_name, quote = FALSE)
      # plot_name <- paste("./Taylor/Bacteria/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
      # par(mar=c(3,5,2,2))
      # p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
      #   geom_point(size=1)  +
      #   geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
      model.lm <- lm(LogVar ~ LogMean, data = all_summary)
      l <- list(beta = format(coef(model.lm)[[2]]),
                r2 = format(summary(model.lm)$r.squared, digits = 3),
                p = format(summary(model.lm)$coefficients[2,4], digits = 3))
      # eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
      # p + theme_bw() +
      #   geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
      #   labs(x="Log mean abundance", y = "Log variance")
      # dev.off()
      # ggsave(plot_name, width = 4, height = 4)
      # dev.new()
      beta_summary = data.frame(Plot=Plot_list[i], Phylum=Phylum_list[phy], Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
      beta_summary_all = rbind(beta_summary_all, beta_summary)
    }
  }
}
write.csv(beta_summary_all, "./Taylor/Bacteria/Phylum/UPARSE_diff_Phylum_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Fungi data
##UPARSE
##Group by taxonomy(Phylum level)
beta_summary_all = data.frame(Plot=character(),	Phylum=character(), Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Fungi/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
  if(file.exists(data_name)){
    subcom = read.csv(data_name, header=T, row.names =1)
    taxonomy = read.csv(paste("./data/Fungi/", Plot_list[i], "_UPARSE_otu_All_taxonomy_abun_freq.csv",sep=""), header=T, row.names =1)
    Phylum_num = as.data.frame(table(taxonomy$Phylum))
    colnames(Phylum_num) = c("Phylum","Number")
    write.csv(Phylum_num, paste("./Taylor/Fungi/Phylum/", Plot_list[i], "_UPARSE_otu_All_taxonomy_Phylum_number.csv",sep=""), quote = F, row.names = F)
    Phylum_list = as.character(Phylum_num[Phylum_num$Number>15,]$Phylum)
    for(phy in 1:length(Phylum_list)){
      subcom_select = subcom[rownames(taxonomy[taxonomy$Phylum==Phylum_list[phy],]),]
      all_summary <- data.frame(Mean = double(), Var = double())
      for(n in 1:nrow(subcom_select)){
        abun <- as.numeric(subcom_select[n,])
        mean_abun <- mean(abun)
        var_abun <- var(abun)
        abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
        rownames(abun_summary) <- rownames(subcom_select)[n]
        all_summary <- rbind(all_summary, abun_summary)
      }
      ##plot real point and fit line
      all_summary = all_summary[all_summary$Variance>0,]
      all_summary["LogMean"] <- log10(all_summary$Mean)
      all_summary["LogVar"] <- log10(all_summary$Variance)
      summary_name <- paste("./Taylor/Fungi/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
      write.csv(all_summary, summary_name, quote = FALSE)
      # plot_name <- paste("./Taylor/Fungi/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
      # par(mar=c(3,5,2,2))
      # p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
      #   geom_point(size=1)  +
      #   geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
      model.lm <- lm(LogVar ~ LogMean, data = all_summary)
      if(is.na(coef(model.lm)[[2]])){
        l <- list(beta = 0, r2 = 0, p = 0)
      }else{
        l <- list(beta = format(coef(model.lm)[[2]]),
                  r2 = format(summary(model.lm)$r.squared, digits = 3),
                  p = format(summary(model.lm)$coefficients[2,4], digits = 3))
      }
      # eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
      # p + theme_bw() +
      #   geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
      #   labs(x="Log mean abundance", y = "Log variance")
      # dev.off()
      # ggsave(plot_name, width = 4, height = 4)
      # dev.new()
      beta_summary = data.frame(Plot=Plot_list[i], Phylum=Phylum_list[phy], Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
      beta_summary_all = rbind(beta_summary_all, beta_summary)
    }
  }
}
write.csv(beta_summary_all, "./Taylor/Fungi/Phylum/UPARSE_diff_Phylum_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Protists data
##UPARSE
##Group by taxonomy(Phylum level)
beta_summary_all = data.frame(Plot=character(),	Phylum=character(), Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  data_name = paste("./data/Protists/", Plot_list[i], "_UPARSE_otu_All_abun.csv",sep="")
  if(file.exists(data_name)){
    subcom = read.csv(data_name, header=T, row.names =1)
    taxonomy = read.csv(paste("./data/Protists/", Plot_list[i], "_UPARSE_otu_All_taxonomy_abun_freq.csv",sep=""), header=T, row.names =1)
    Phylum_num = as.data.frame(table(taxonomy$Phylum))
    colnames(Phylum_num) = c("Phylum","Number")
    write.csv(Phylum_num, paste("./Taylor/Protists/Phylum/", Plot_list[i], "_UPARSE_otu_All_taxonomy_Phylum_number.csv",sep=""), quote = F, row.names = F)
    Phylum_list = as.character(Phylum_num[Phylum_num$Number>15,]$Phylum)
    for(phy in 1:length(Phylum_list)){
      subcom_select = subcom[rownames(taxonomy[taxonomy$Phylum==Phylum_list[phy],]),]
      all_summary <- data.frame(Mean = double(), Var = double())
      for(n in 1:nrow(subcom_select)){
        abun <- as.numeric(subcom_select[n,])
        mean_abun <- mean(abun)
        var_abun <- var(abun)
        abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
        rownames(abun_summary) <- rownames(subcom_select)[n]
        all_summary <- rbind(all_summary, abun_summary)
      }
      ##plot real point and fit line
      all_summary = all_summary[all_summary$Variance>0,]
      all_summary["LogMean"] <- log10(all_summary$Mean)
      all_summary["LogVar"] <- log10(all_summary$Variance)
      summary_name <- paste("./Taylor/Protists/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
      write.csv(all_summary, summary_name, quote = FALSE)
      # plot_name <- paste("./Taylor/Protists/Phylum/",Plot_list[i],"_", Phylum_list[phy],"_UPARSE_all_otus_mean_var_point_fit_line.pdf",sep="")
      # par(mar=c(3,5,2,2))
      # p <- ggplot(data = all_summary, aes(x = LogMean, y = LogVar)) +
      #   geom_point(size=1)  +
      #   geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",linewidth=1) 
      model.lm <- lm(LogVar ~ LogMean, data = all_summary)
      if(is.na(coef(model.lm)[[2]])){
        l <- list(beta = 0, r2 = 0, p = 0)
      }else{
        l <- list(beta = format(coef(model.lm)[[2]]),
                  r2 = format(summary(model.lm)$r.squared, digits = 3),
                  p = format(summary(model.lm)$coefficients[2,4], digits = 3))
      }
      # eq <- substitute(italic(B)~"="~beta~","~italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
      # p + theme_bw() +
      #   geom_text(aes(x = 0.5*(min(LogMean)+max(LogMean)), y = 1.1*max(LogVar), label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
      #   labs(x="Log mean abundance", y = "Log variance")
      # dev.off()
      # ggsave(plot_name, width = 4, height = 4)
      # dev.new()
      beta_summary = data.frame(Plot=Plot_list[i], Phylum=Phylum_list[phy], Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
      beta_summary_all = rbind(beta_summary_all, beta_summary)
    }
  }
}
write.csv(beta_summary_all, "./Taylor/Protists/Phylum/UPARSE_diff_Phylum_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)


##Nematoda data
##UPARSE
##Group by taxonomy(Trophic level)
Trophic_list = c("b","f","h","o","p")
info = lazyopen("./data/summary of environmental data2010-2020.txt")
Plot_list = unique(info$plot)
beta_summary_all = data.frame(Plot=character(),	Trophic=character(), Intercept = double(), Beta=double(), R2=double(), P=double())
for(i in 1:length(Plot_list)){
  for(j in 1:length(Trophic_list)){
    data_name = paste("./data/Nematoda/Trophic/Trophic_",Trophic_list[j],"_",Plot_list[i], "_UPARSE_otu_All.csv",sep="")
    all_summary_all = data.frame(Mean = double(), Variance = double(), LogMean = double(), LogVar=double(), Trophic = character())
    if(file.exists(data_name)){
      subcom = read.csv(data_name, header=T, row.names =1)
      if(nrow(subcom)>5){
        all_summary <- data.frame(Mean = double(), Var = double())
        for(n in 1:nrow(subcom)){
          abun <- as.numeric(subcom[n,])
          mean_abun <- mean(abun)
          var_abun <- var(abun)
          abun_summary <- data.frame(Mean = mean_abun, Variance = var_abun)
          rownames(abun_summary) <- rownames(subcom)[n]
          all_summary <- rbind(all_summary, abun_summary)
        }
        ##plot real point and fit line
        all_summary = all_summary[all_summary$Variance>0,]
        all_summary["LogMean"] <- log10(all_summary$Mean)
        all_summary["LogVar"] <- log10(all_summary$Variance)
        all_summary["Trophic"] <- Trophic_list[j]
        all_summary_all = rbind(all_summary_all, all_summary)
        #summary_name <- paste("./Taylor/Nematoda/Trophic/",Plot_list[i],"_", Trophic_list[phy],"_UPARSE_all_otus_mean_var_summary.csv", sep='')
        #write.csv(all_summary, summary_name, quote = FALSE)
        model.lm <- lm(LogVar ~ LogMean, data = all_summary)
        if(is.na(coef(model.lm)[[1]])|is.na(coef(model.lm)[[2]])){
          l <- list(intercept = 0, beta = 0, r2 = 0, p = 1)
        }else{
          l <- list(intercept = format(coef(model.lm)[[1]], digits = 3),
                    beta = format(coef(model.lm)[[2]]),
                    r2 = format(summary(model.lm)$r.squared, digits = 3),
                    p = format(summary(model.lm)$coefficients[2,4], digits = 3))
        }
        beta_summary = data.frame(Plot=Plot_list[i], Trophic=Trophic_list[j], Intercept= as.numeric(l$intercept), Beta=as.numeric(l$beta), R2=as.numeric(l$r2), P=as.numeric(l$p))
        beta_summary_all = rbind(beta_summary_all, beta_summary)
      }
    }
    write.csv(all_summary_all, paste("./Taylor/Nematoda/Trophic/Trophic_",Trophic_list[j],"_",Plot_list[i],"_UPARSE_all_otus_mean_var_summary.csv", sep=''), quote = FALSE)
  }
}
write.csv(beta_summary_all, "./Taylor/Nematoda/Trophic/UPARSE_diff_Trophic_all_otus_Taylor_beta.csv", quote = FALSE, row.names = FALSE)
