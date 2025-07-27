rm(list=ls())
library(iCAMP)
packageVersion("iCAMP")
library(ieggr)
library(castor)
library(ape)
library(picante)

##calculate MPDik of each OTU in each sample 
mpdi.big <- function (comm, pd.desc = "pd.desc", pd.spname = NULL, pd.wd = getwd(), 
                      sp.limit = 10000, spname.check = TRUE, abundance.weighted = TRUE) 
{
  requireNamespace("bigmemory")
  if (spname.check) {
    if (is.null(pd.spname)) 
      pd.spname = lazyopen(file = paste0(pd.wd, "/pd.taxon.name.csv"))[, 
                                                                       1]
    check.sp = match.name(name.check = pd.spname, cn.list = list(comm = comm))
    comm = check.sp$comm
  }
  pdbig.id = match(colnames(comm), pd.spname)
  if (!abundance.weighted) {
    comm[comm > 0] = 1
    num = rowSums(comm)
  }
  comm = comm/rowSums(comm)
  comm = as.matrix(comm)
  pd = try(bigmemory::attach.big.matrix(dget(paste0(pd.wd, 
                                                    "/", pd.desc))))
  if (inherits(pd, "try-error")) {
    pd = bigmemory::attach.big.matrix(paste0(pd.wd, "/", 
                                             pd.desc))
  }
  if (nrow(pd) <= sp.limit) {
    comd = comm %*% pd[pdbig.id, pdbig.id]
  }else {
    N.int = floor(sp.limit * sp.limit/ncol(comm))
    N.num = ceiling(ncol(comm)/N.int)
    N.ser = ((1:N.num) - 1) * N.int + 1
    N.ser = c(N.ser, ncol(pd) + 1)
    comd = comm
    for (i in 1:N.num) {
      message("Now begin i=", i, " in ", N.num, ", with memory of ", 
              memory.size(), " Mb. ", date())
      comd[, N.ser[i]:(N.ser[i + 1] - 1)] = comm %*% pd[pdbig.id, 
                                                        pdbig.id[N.ser[i]:(N.ser[i + 1] - 1)]]
      message("temp memory is ", memory.size(), " Mb.")
      gc()
    }
  }
  gc()
  res = comd * comm
  colnames(res) = colnames(comm)
  res
}
##calculate MNTDik of each OTU in each sample 
mntdi.big <- function (comm, pd.desc = "pd.desc", pd.spname = NULL, pd.wd = getwd(), 
                       spname.check = FALSE, abundance.weighted = TRUE, memory.G = 50) 
{
  requireNamespace("bigmemory")
  if (.Platform$OS.type == "windows") {
    if (utils::memory.limit() < memory.G * 1024) {
      memotry = try(utils::memory.limit(size = memory.G * 
                                          1024), silent = TRUE)
      if (inherits(memotry, "try-error")) {
        warning(memotry[1])
      }
    }
  }
  if (spname.check) {
    if (is.null(pd.spname)) 
      pd.spname = lazyopen(file = paste0(pd.wd, "/pd.taxon.name.csv"))[, 
                                                                       1]
    check.sp = match.name(name.check = pd.spname, cn.list = list(comm = comm))
    comm = check.sp$comm
  }
  pdbig.id = match(colnames(comm), pd.spname)
  pd = try(bigmemory::attach.big.matrix(dget(paste0(pd.wd, 
                                                    "/", pd.desc))))
  if (inherits(pd, "try-error")) {
    pd = bigmemory::attach.big.matrix(paste0(pd.wd, "/", 
                                             pd.desc))
  }
  N = nrow(comm)
  gc()
  if (!abundance.weighted) {
    comm[comm > 0] = 1
    num = rowSums(comm)
  }
  min.d = matrix(0, nrow = N, ncol = ncol(comm))
  for (i in 1:N) {
    id = which(comm[i, ] > 0)
    pdx = pd[pdbig.id[id], pdbig.id[id]]
    diag(pdx) = NA
    min.d[i, id] = apply(pdx, 2, min, na.rm = TRUE)
    gc()
  }
  comm.p = comm/rowSums(comm)
  res = min.d * comm.p
  res
}

##Archaea
if(!file.exists("Archaea/")){
  dir.create("Archaea")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample7842_otutab_Archaea_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_16S_archaea_2010to2020_resample7842.nwk")
pdist.out = iCAMP::pdist.big(mytree, wd=paste0(getwd(),"/Archaea"))
save(pdist.out, file = "./Archaea/16Sotu_2010to2020_pdist_out.rda")
mycomm_mpd = iCAMP::mpd.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                            pd.wd = pdist.out$pd.wd, spname.check=TRUE)
write.csv(as.data.frame(mycomm_mpd), "./Archaea/MPD_in_each_sample_16Sotu_2010to2020.csv", quote = F)
each_otu_mpd = mpdi.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                        pd.wd = pdist.out$pd.wd, spname.check=TRUE)
otu_MPDik = each_otu_mpd/mycomm_mpd
write.csv(otu_MPDik, "./Archaea/MPD_weighted_16Sotu_2010to2020.csv", quote = F)


##Bacteria
if(!file.exists("Bacteria/")){
  dir.create("Bacteria")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample26890_otutab_Bacteria_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Bacteria_2010to2020_resample26890.nwk")
pdist.out = iCAMP::pdist.big(mytree, wd=paste0(getwd(),"./Bacteria"))
save(pdist.out, file = "./Bacteria/16Sotu_2010to2020_pdist_out.rda")
mycomm_mpd = iCAMP::mpd.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                            pd.wd = pdist.out$pd.wd, spname.check=TRUE)
write.csv(as.data.frame(mycomm_mpd), "./Bacteria/MPD_in_each_sample_16Sotu_2010to2020withreseq0608.csv", quote = F)
each_otu_mpd = mpdi.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                        pd.wd = pdist.out$pd.wd, spname.check=TRUE)
otu_MPDik = each_otu_mpd/mycomm_mpd
write.csv(otu_MPDik, "./Bacteria/MPD_weighted_16Sotu_2010to2020withreseq0608.csv", quote = F)


##Fungi
if(!file.exists("Fungi/")){
  dir.create("Fungi")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample541_otutab_Fungi_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Fungi_2010to2020_resample541.nwk")
pdist.out = iCAMP::pdist.big(mytree, wd=paste0(getwd(),"/Fungi"))
save(pdist.out, file = "./Fungi/18Sotu_2010to2020_pdist_out.rda")
mycomm_mpd = iCAMP::mpd.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                            pd.wd = pdist.out$pd.wd, spname.check=TRUE)
write.csv(as.data.frame(mycomm_mpd), "./Fungi/MPD_in_each_sample_18Sotu_2010to2020withreseq0608.csv", quote = F)
each_otu_mpd = mpdi.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                        pd.wd = pdist.out$pd.wd, spname.check=TRUE)
otu_MPDik = each_otu_mpd/mycomm_mpd
write.csv(otu_MPDik, "./Fungi/MPD_weighted_18Sotu_2010to2020withreseq0608.csv", quote = F)


##Protists
if(!file.exists("Protists/")){
  dir.create("Protists")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample1064_otutab_Protists_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Protists_2010to2020_resample1064.nwk")
match.phylo.tree = match.phylo.data(mytree, mycomm)
mycomm = match.phylo.tree$data
mytree = match.phylo.tree$phy
pdist.out = iCAMP::pdist.big(mytree, wd=paste0(getwd(),"/Protists"))
save(pdist.out, file = "./Protists/18Sotu_2010to2020_pdist_out.rda")
mycomm_mpd = iCAMP::mpd.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                            pd.wd = pdist.out$pd.wd, spname.check=TRUE)
write.csv(as.data.frame(mycomm_mpd), "./Protists/MPD_in_each_sample_18Sotu_2010to2020withreseq0608.csv", quote = F)
each_otu_mpd = mpdi.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                        pd.wd = pdist.out$pd.wd, spname.check=TRUE)
otu_MPDik = each_otu_mpd/mycomm_mpd
write.csv(otu_MPDik, "./Protists/MPD_weighted_18Sotu_2010to2020withreseq0608.csv", quote = F)


##Nematoda
if(!file.exists("Nematoda/")){
  dir.create("Nematoda")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample296_otutab_Nematoda_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Nematoda_2010to2020_resample296.nwk")
pdist.out = iCAMP::pdist.big(mytree, wd=paste0(getwd(),"/Nematoda"))
save(pdist.out, file = "./Nematoda/18Sotu_2010to2020_pdist_out.rda")
mycomm_mpd = iCAMP::mpd.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                            pd.wd = pdist.out$pd.wd, spname.check=TRUE)
write.csv(as.data.frame(mycomm_mpd), "./Nematoda/MPD_in_each_sample_18Sotu_2010to2020withreseq0608.csv", quote = F)
each_otu_mpd = mpdi.big(t(mycomm), pd.desc = pdist.out$pd.file, pd.spname=pdist.out$tip.label,
                        pd.wd = pdist.out$pd.wd, spname.check=TRUE)
otu_MPDik = each_otu_mpd/mycomm_mpd
write.csv(otu_MPDik, "./Nematoda/MPD_weighted_18Sotu_2010to2020withreseq0608.csv", quote = F)
