rm(list=ls())
require(ggplot2)
require(ieggr)
library(dplyr)
require(reshape2)
library(scatterplot3d)
library(hrbrthemes)
##Taxonomic
##UPARSE
##Archaea data
Archaea_beta_compare = read.csv("./Taylor/Archaea/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Archaea_beta_compare$Taxa = "Archaea"
##Bacteria data
Bacteria_beta_compare = read.csv("./Taylor/Bacteria/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Bacteria_beta_compare$Taxa = "Bacteria"
##Fungi data
Fungi_beta_compare = read.csv("./Taylor/Fungi/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Fungi_beta_compare$Taxa = "Fungi"
##Protists data
Protists_beta_compare = read.csv("./Taylor/Protists/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Protists_beta_compare$Taxa = "Protists"
##Nematoda data
Nematoda_beta_compare = read.csv("./Taylor/Nematoda/compare/UPARSE_all_otus_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Nematoda_beta_compare$Taxa = "Nematoda"
Taxonomic_beta_compare = rbind(Archaea_beta_compare, Bacteria_beta_compare, Fungi_beta_compare, Protists_beta_compare, Nematoda_beta_compare)
Taxonomic_beta_compare$Group = "tTTPL exponent"
##Phylogenetic
##MPD
##Archaea data
Archaea_beta_compare = read.csv(".Taylor/Archaea/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Archaea_beta_compare$Taxa = "Archaea"
##Bacteria data
Bacteria_beta_compare = read.csv(".Taylor/Bacteria/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Bacteria_beta_compare$Taxa = "Bacteria"
##Fungi data
Fungi_beta_compare = read.csv(".Taylor/Fungi/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Fungi_beta_compare$Taxa = "Fungi"
##Protists data
Protists_beta_compare = read.csv(".Taylor/Protists/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Protists_beta_compare$Taxa = "Protists"
##Nematoda data
Nematoda_beta_compare = read.csv(".Taylor/Nematoda/compare/MPD_weighted_otu_Taylor_beta_Treatments_control_compare.csv", header=T, row.names = 1)
Nematoda_beta_compare$Taxa = "Nematoda"
Phylogenetic_beta_compare = rbind(Archaea_beta_compare, Bacteria_beta_compare, Fungi_beta_compare, Protists_beta_compare, Nematoda_beta_compare)
Phylogenetic_beta_compare$Group = "pTTPL exponent"
all_beta_compare = rbind(Taxonomic_beta_compare, Phylogenetic_beta_compare)
write.csv(all_beta_compare, "./Figures/Figure2a_data.csv", quote = F)
##single factors
all_beta_compare = lazyopen("./Figures/Figure2a_data.csv")
beta_compare_select = all_beta_compare[all_beta_compare$Treatment%in%c("W","H","D","C"),]
beta_compare_select$Treatment = factor(beta_compare_select$Treatment, levels = c("W","H","D","C"))
beta_compare_select$Taxa = factor(beta_compare_select$Taxa, levels = c("Bacteria","Archaea","Fungi","Protists","Nematoda"))
beta_compare_select$Group = factor(beta_compare_select$Group, levels = c("tTTPL exponent","pTTPL exponent"))
par(mar=c(3,5,2,2))
p<- ggplot(beta_compare_select, aes(x=Treatment, y=Cohen.d, group=Taxa, fill=Taxa)) + 
  geom_col(position=position_dodge(width=0.9), alpha=0.6) +
  scale_fill_manual(values = c("#E69F00", "#8B9E79", "#009E73", "#D55E00","#0072B9"), limits = c("Bacteria","Archaea","Fungi","Protists","Nematoda"))+
  ylim(-4,4)
p + labs(x=NULL, y = "Cohen's d") +
  theme_bw() +
  theme(axis.text = element_text(face = "plain", size = 12, color = "black"),
        axis.title = element_text(face = "plain", size = 12, color = "black"),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'))+
  #theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  theme(legend.position = "right")+
  facet_wrap(Group~., ncol = 2, scales="free")
dev.off()
ggsave("./Figures/Figure2a_plot.pdf", width = 9, height = 3)

###3D plot
##combined factors
##tTTPL
all_beta = lazyopen("./Figures/Figure2a_data.csv")
tTTPL_beta = all_beta[all_beta$Group = "tTTPL exponent",]
tTTPL_beta = tTTPL_beta[,c("Treatment","Taxa","group","colorgroup","Cohen.d","RR.P.FillNA","Sig.RR.FillNA")]
tTTPL_beta$type = "tTTPL"
write.csv(tTTPL_beta, "Figure2b_data.csv", quote = F)
all_beta = lazyopen("Figure2b_data.csv")
select_beta = all_beta[all_beta$Taxa %in% c("Bacteria","Protists","Nematoda"),]
select_tTTPL = dcast(select_beta, Treatment ~ Taxa, value.var="Cohen.d")
select_tTTPL$colorgroup = 1
select_tTTPL[select_tTTPL$Treatment%in%c("W","WH","WD","WC","WHC","WDC"),]$colorgroup = 2
##plot
par(mar = c(5.1, 4.1, 4.1, 7.1))
colors0 <-  c("black","#E9212C")
colors <- colors0[select_tTTPL$colorgroup]
p3d <- scatterplot3d(select_tTTPL[,2:4], pch = 16, color=colors, angle=40, grid=F,
                     xlab = "Bacteria d", ylab = "Protists d",zlab = "Nematoda d")
legend("right", legend = c("unwarm", "warm"), col =  c("black","#E9212C"), pch = 16)
text(p3d$xyz.convert(select_tTTPL[, 2:4]), labels = select_tTTPL$Treatment,
     cex= 1, col = "black")
##Archaea~Fungi
select_beta = all_beta[all_beta$Taxa %in% c("Archaea","Fungi"),]
select_tTTPL = dcast(select_beta, Treatment ~ Taxa, value.var="Cohen.d")
select_tTTPL$colorgroup = "unwarm"
select_tTTPL[select_tTTPL$Treatment%in%c("W","WH","WD","WC","WHC","WDC"),]$colorgroup = "warm"
par(mar=c(3,5,2,2))
p<- ggplot(select_tTTPL, aes(x=Fungi, y=Archaea, color=colorgroup)) + 
  geom_point(size=2) +
  geom_text(aes(label = Treatment), vjust = -1, hjust = 0.5, size=3) +
  scale_color_manual(values = c("black", "red"), limits = c("unwarm","warm"))+
  scale_x_continuous(breaks = seq(-0.5, 3, by = 0.5)) +
  scale_y_continuous(breaks = seq(-0.5, 1.5, by = 0.5)) +
  geom_vline(aes(xintercept=0), colour = "black", lty="dashed", linewidth =0.5) +
  geom_hline(aes(yintercept=0), colour = "black", lty="dashed", linewidth =0.5)
p + labs(x="Fungi d", y = "Archaea d") +
  theme_bw() +
  theme(axis.text = element_text(face = "plain", size = 12, color = "black"),
        axis.title = element_text(face = "plain", size = 12, color = "black"),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position = "none")
dev.off()
ggsave("./FigureS5a.pdf", width = 4, height = 4)


##pTTPL
all_beta = lazyopen("./Figures/Figure2a_data.csv")
pTTPL_beta = all_beta[all_beta$Group = "pTTPL exponent",]
pTTPL_beta = pTTPL_beta[,c("Treatment","Taxa","group","colorgroup","Cohen.d","RR.P.FillNA","Sig.RR.FillNA")]
pTTPL_beta$type = "pTTPL"
write.csv(pTTPL_beta, "Figure2c_data.csv", quote = F)
all_beta = lazyopen("Figure2c_data.csv")
select_beta = all_beta[all_beta$Taxa %in% c("Bacteria","Protists","Nematoda"),]
select_pTTPL = dcast(select_beta, Treatment ~ Taxa, value.var="Cohen.d")
select_pTTPL$colorgroup = 1
select_pTTPL[select_pTTPL$Treatment%in%c("W","WH","WD","WC","WHC","WDC"),]$colorgroup = 2
##plot
par(mar = c(5.1, 4.1, 4.1, 7.1))
colors0 <-  c("black","#E9212C")
colors <- colors0[select_pTTPL$colorgroup]
p3d <- scatterplot3d(select_pTTPL[,2:4], pch = 16, color=colors, angle=40, grid=F,
                     xlab = "Bacteria d", ylab = "Protists d",zlab = "Nematoda d")
legend("right", legend = c("unwarm", "warm"), col =  c("black","#E9212C"), pch = 16)
text(p3d$xyz.convert(select_pTTPL[, 2:4]), labels = select_pTTPL$Treatment,
     cex= 1, col = "black")
##Archaea~Fungi
select_beta = all_beta[all_beta$Taxa %in% c("Archaea","Fungi"),]
select_pTTPL = dcast(select_beta, Treatment ~ Taxa, value.var="Cohen.d")
select_pTTPL$colorgroup = "unwarm"
select_pTTPL[select_pTTPL$Treatment%in%c("W","WH","WD","WC","WHC","WDC"),]$colorgroup = "warm"
par(mar=c(3,5,2,2))
p<- ggplot(select_pTTPL, aes(x=Fungi, y=Archaea, color=colorgroup)) + 
  geom_point(size=2) +
  geom_text(aes(label = Treatment), vjust = -1, hjust = 0.5, size=3) +
  scale_color_manual(values = c("black", "red"), limits = c("unwarm","warm"))+
  scale_x_continuous(breaks = seq(-0.5, 3, by = 0.5)) +
  scale_y_continuous(breaks = seq(-0.5, 2, by = 0.5)) +
  geom_vline(aes(xintercept=0), colour = "black", lty="dashed", linewidth =0.5) +
  geom_hline(aes(yintercept=0), colour = "black", lty="dashed", linewidth =0.5)
p + labs(x="Fungi d", y = "Archaea d") +
  theme_bw() +
  theme(axis.text = element_text(face = "plain", size = 12, color = "black"),
        axis.title = element_text(face = "plain", size = 12, color = "black"),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.position = "none")
dev.off()
ggsave("./FigureS5b.pdf", width = 4, height = 4)
