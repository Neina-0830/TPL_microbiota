# Phylogenetic distance versus differences in group-level beta and warming response.
source(file.path(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),"code/analysis_helpers.R"))
setwd(project_root());require_packages("ape")
meta <- read.csv("data/prepared/environment_metadata.csv",check.names=FALSE)
spec <- data.frame(Group=c("Bacteria","Archaea","Fungi","Protists","Nematoda"),
  Quantile=c(.01,.4,.2,.1,.4),Focal=c(1000,Inf,Inf,500,Inf),Seed=20261021:20261025)
N_PERM <- as.integer(Sys.getenv("TPL_PHYLOGROUP_N_PERM",unset="9999"))
out <- "results/06_phylogroups";dir.create(out,recursive=TRUE,showWarnings=FALSE)
summary <- list()
for(i in seq_len(nrow(spec))) {
  g <- spec$Group[i];set.seed(spec$Seed[i])
  x <- as.matrix(read.csv(file.path("data/prepared/full_relative_abundance",paste0(g,"_relative_abundance_2010to2020.csv")),row.names=1,check.names=FALSE))
  # Relative abundance has already been calculated from the complete community.
  treefile <- list.files("data",pattern=paste0("^rooted_tree_otus_",g,"_"),full.names=TRUE)
  tree <- ape::read.tree(treefile)
  taxa <- intersect(rownames(x)[rowMeans(x>0)>=.10],tree$tip.label)
  tree <- ape::keep.tip(tree,taxa);x <- x[tree$tip.label,,drop=FALSE]
  D <- ape::cophenetic.phylo(tree)
  focal <- if(is.finite(spec$Focal[i]) && length(taxa)>spec$Focal[i]) sample(rownames(x),spec$Focal[i]) else rownames(x)
  fd <- D[focal,,drop=FALSE]
  # Exclude focal-to-self distances explicitly before selecting the threshold.
  distances <- unlist(lapply(focal,function(f) D[f,setdiff(colnames(D),f)]),use.names=FALSE)
  threshold <- quantile(distances,spec$Quantile[i],names=FALSE)
  plots <- unique(meta$plot)
  beta <- matrix(NA_real_,length(focal),length(plots),dimnames=list(focal,plots))
  membership <- lapply(focal,function(f) colnames(D)[D[f,]<=threshold])
  for(p in plots) {
    md <- meta[meta$plot==p,,drop=FALSE];md <- md[order(md$year),]
    a <- x[,md$Sample_name,drop=FALSE];mu <- rowMeans(a);va <- apply(a,1,var)
    for(j in seq_along(focal)) {
      taxa_j <- membership[[j]];valid <- taxa_j[mu[taxa_j]>0 & va[taxa_j]>0 & is.finite(mu[taxa_j]) & is.finite(va[taxa_j])]
      if(ncol(a)>=8 && length(valid)>=20) beta[j,p] <- tpl_fit(mu[valid],va[valid])["Beta"]
    }
  }
  beta <- beta[apply(beta,1,function(z) all(is.finite(z))),,drop=FALSE]
  if(nrow(beta)<3) stop("Too few complete phylogroups: ",g)
  focal <- rownames(beta);df <- D[focal,focal,drop=FALSE]
  values <- list(Mean_beta=rowMeans(beta))
  block_deltas <- lapply(unique(meta$block),function(b) {
    md <- unique(meta[meta$block==b,c("plot","warm")]);
    rowMeans(beta[,md$plot[md$warm==1],drop=FALSE])-rowMeans(beta[,md$plot[md$warm==0],drop=FALSE])
  })
  values$Warming_response <- rowMeans(do.call(cbind,block_deltas))
  write_table(data.frame(Focal=focal,beta,check.names=FALSE),file.path(out,paste0(g,"_subplot_betas.csv")))
  tri <- lower.tri(df);rd <- rank(df[tri]);rd <- rd-mean(rd);nd <- sqrt(sum(rd^2))
  for(analysis in names(values)) {
    v <- values[[analysis]];dv <- abs(outer(v,v,"-"));rv <- rank(dv[tri]);rv <- rv-mean(rv)
    rankmat <- matrix(0,length(v),length(v));rankmat[tri] <- rv;rankmat <- rankmat+t(rankmat)
    denom <- nd*sqrt(sum(rv^2));observed <- sum(rd*rv)/denom
    null <- replicate(N_PERM,{p <- sample.int(length(v));sum(rd*rankmat[p,p][tri])/denom})
    summary[[paste(g,analysis)]] <- data.frame(Group=g,Analysis=analysis,N_focal=length(v),
      Threshold=threshold,Spearman_rho=observed,P=(1+sum(null>=observed))/(1+N_PERM))
  }
}
s <- do.call(rbind,summary);s$P_Holm <- ave(s$P,s$Analysis,FUN=function(p) p.adjust(p,"holm"))
write_table(s,file.path(out,"phylogroup_tests.csv"))
