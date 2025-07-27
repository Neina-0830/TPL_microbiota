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

generate_null_mpdik_tipshuffle <- function(comm, tree, pdist.out, nrep = 999, output_dir = NULL) {
  common_otus <- intersect(rownames(comm), pdist.out$tip.label)
  comm <- comm[common_otus, , drop = FALSE]
  tree <- keep.tip(tree, common_otus)
  for (i in 1:nrep) {
    cat("Null model iteration", i, "of", nrep, "\n")
    # shuffle tip labels
    shuffled_tree <- tree
    shuffled_tree$tip.label <- sample(tree$tip.label)
    # new temp pd file
    temp_subdir <- file.path(tempdir(), paste0("null_shuffle_", i))
    dir.create(temp_subdir, showWarnings = FALSE)
    shuffled_pdist <- iCAMP::pdist.big(shuffled_tree, wd = temp_subdir)
    comm_mpd = iCAMP::mpd.big(t(comm), pd.desc = shuffled_pdist$pd.file, pd.spname=shuffled_pdist$tip.label,
                              pd.wd = shuffled_pdist$pd.wd, spname.check=TRUE)
    # calculate MPDik
    each_otu_mpd <- mpdi.big(t(comm), pd.desc = shuffled_pdist$pd.file, pd.spname = shuffled_pdist$tip.label,
                             pd.wd = shuffled_pdist$pd.wd, spname.check = TRUE)
    otu_MPDik = each_otu_mpd/comm_mpd
    if (!is.null(output_dir)) {
      write.csv(otu_MPDik,
                file = file.path(output_dir, paste0("MPDik_null_", i, ".csv")),
                quote = FALSE)
    }
    dirs_to_remove <- list.dirs(tempdir(), full.names = TRUE, recursive = FALSE)
    dirs_to_remove <- dirs_to_remove[grepl("null_shuffle_", basename(dirs_to_remove))]
    unlink(dirs_to_remove, recursive = TRUE)
  }
}

##Archaea
if(!file.exists("Archaea/Null")){
  dir.create("Archaea/Null")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample7842_otutab_Archaea_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_16S_archaea_2010to2020_resample7842.nwk")
load("./Archaea/16Sotu_2010to2020_pdist_out.rda")
my_dir = "Archaea/Null"
generate_null_mpdik_tipshuffle(mycomm, mytree, pdist.out, nrep = 100, output_dir=my_dir)


##Bacteria
if(!file.exists("Bacteria/Null")){
  dir.create("Bacteria/Null")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample26890_otutab_Bacteria_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Bacteria_2010to2020_resample26890.nwk")
load("./Bacteria/16Sotu_2010to2020_pdist_out.rda")
my_dir = "Bacteria/Null"
generate_null_mpdik_tipshuffle(mycomm, mytree, pdist.out, nrep = 100, output_dir=my_dir)


##Fungi
if(!file.exists("Fungi/Null")){
  dir.create("Fungi/Null")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample541_otutab_Fungi_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Fungi_2010to2020_resample541.nwk")
load("./Fungi/18Sotu_2010to2020_pdist_out.rda")
my_dir = "Fungi/Null"
generate_null_mpdik_tipshuffle(mycomm, mytree, pdist.out, nrep = 100, output_dir=my_dir)


##Protists
if(!file.exists("Protists/Null")){
  dir.create("Protists/Null")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./rawdata/resample1064_otutab_Protists_2010to2020.csv")
mytree = read.tree("./rawdata/rooted_tree_otus_Protists_2010to2020_resample1064.nwk")
load("./Protists/18Sotu_2010to2020_pdist_out.rda")
my_dir = "Protists/Null"
generate_null_mpdik_tipshuffle(mycomm, mytree, pdist.out, nrep = 100, output_dir=my_dir)


##Nematoda
if(!file.exists("Nematoda/Null")){
  dir.create("Nematoda/Null")
}
##calculate mean pairwise distance of all samples
mycomm = lazyopen("./data/resample296_otutab_Nematoda_2010to2020.csv")
mytree = read.tree("./data/rooted_tree_otus_Nematoda_2010to2020_resample296.nwk")
load("./Nematoda/18Sotu_2010to2020_pdist_out.rda")
my_dir = "Nematoda/Null"
generate_null_mpdik_tipshuffle(mycomm, mytree, pdist.out, nrep = 100, output_dir=my_dir)
