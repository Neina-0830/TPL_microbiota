# Relative effects of single-process removal on the TPL exponent.
source(file.path(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),"code/analysis_helpers.R"))
setwd(project_root())
ROOT <- Sys.getenv("TPL_SGLV_WORK_DIR",unset="results/17_sGLV")
INPUT <- Sys.getenv("TPL_SGLV_BETA_FILE",unset=file.path(ROOT,"Taylor/process_removal_top10/TPL_estimates.csv"))
OUT <- Sys.getenv("TPL_SGLV_IMPORTANCE_OUT",unset=file.path(ROOT,"process_importance"))
PLOT_ORDER <- c("immigration","noise","growth","interaction")
ALPHA <- .05;ZERO_TOL <- 1e-12
BOOT_REPS <- as.integer(Sys.getenv("TPL_BOOT_REPS",unset="20000"))
set.seed(20261017)
compact_letters <- function(pmat,means,alpha=ALPHA) {
  nm <- names(means); k <- length(nm)
  if(!all(dim(pmat)==c(k,k))) stop("Invalid P matrix.")
  if(any(!is.finite(pmat))) return(setNames(rep("",k),nm))
  sig <- pmat<alpha; diag(sig) <- FALSE
  cliques <- list()
  for(mask in seq_len(2^k-1L)) {
    ids <- which(as.logical(intToBits(mask)[seq_len(k)]))
    if(!any(sig[ids,ids,drop=FALSE])) cliques[[length(cliques)+1L]] <- ids
  }
  keep <- vapply(seq_along(cliques),function(i)
    !any(vapply(seq_along(cliques),function(j)
      length(cliques[[j]])>length(cliques[[i]]) &&
      all(cliques[[i]] %in% cliques[[j]]),logical(1))),logical(1))
  cliques <- cliques[keep]
  order_process <- order(-means,seq_along(means))
  process_rank <- match(seq_len(k),order_process)
  oo <- order(vapply(cliques,function(ii) min(process_rank[ii]),numeric(1)),
              -vapply(cliques,length,integer(1)))
  cliques <- cliques[oo]
  lab <- setNames(rep("",k),nm)
  for(i in seq_along(cliques)) lab[cliques[[i]]] <- paste0(lab[cliques[[i]]],letters[i])
  for(i in seq_len(k-1L)) for(j in (i+1L):k) {
    shared <- length(intersect(strsplit(lab[i],"")[[1]],strsplit(lab[j],"")[[1]]))>0
    if(shared==sig[i,j]) stop("Compact-letter consistency check failed.")
  }
  lab
}

dunn_pairs <- function(values,groups,levels=PLOT_ORDER) {
  groups <- as.character(groups); N <- length(values)
  ranks <- rank(values,ties.method="average")
  nr <- table(factor(groups,levels=levels))
  mr <- tapply(ranks,factor(groups,levels=levels),mean)
  tt <- as.numeric(table(values))
  vv <- N*(N+1)/12-sum(tt^3-tt)/(12*(N-1))
  pairs <- combn(levels,2,simplify=FALSE)
  do.call(rbind,lapply(pairs,function(pp) {
    z <- if(vv>ZERO_TOL) unname((mr[pp[1]]-mr[pp[2]])/
                       sqrt(vv*(1/nr[pp[1]]+1/nr[pp[2]]))) else 0
    data.frame(Process1=pp[1],Process2=pp[2],Statistic=z,P_raw=2*pnorm(-abs(z)))
  }))
}


x <- read.csv(INPUT,check.names=FALSE)
keys <- c("Group","TPL","Condition","Rank","Seed","Scenario")
if(anyDuplicated(x[keys]) || any(!is.finite(x$Beta))) stop("Duplicate or non-finite simulation estimates.")
expected <- c("full",paste0(c("immigration","growth","interaction","noise"),"_removed"))
for(panel in split(x,interaction(x$Group,x$TPL,x$Condition,drop=TRUE))) {
  grid <- expand.grid(Rank=1:10,Seed=0:2,Scenario=expected,stringsAsFactors=FALSE)
  if(nrow(panel)!=nrow(grid) || !setequal(panel$Rank,1:10) || !setequal(panel$Seed,0:2) || !setequal(panel$Scenario,expected)) stop("Incomplete full/removal simulations.")
}
av <- aggregate(Beta~Group+TPL+Condition+Rank+Scenario,x,mean)
full <- av[av$Scenario=="full",c("Group","TPL","Condition","Rank","Beta")];names(full)[5] <- "Beta_full"
removed <- av[av$Scenario!="full",];names(removed)[names(removed)=="Beta"] <- "Beta_removed"
effects <- merge(removed,full,by=c("Group","TPL","Condition","Rank"))
effects$Process <- sub("_removed$","",effects$Scenario)
effects$Effect <- effects$Beta_full-effects$Beta_removed
# Normalize separately within every condition and parameter set.
effects$Denominator <- ave(abs(effects$Effect),interaction(effects$Group,effects$TPL,effects$Condition,effects$Rank),FUN=sum)
effects$I_k <- ifelse(effects$Denominator>ZERO_TOL,abs(effects$Effect)/effects$Denominator,NA_real_)
effects$Signed_I_k <- ifelse(effects$Denominator>ZERO_TOL,effects$Effect/effects$Denominator,NA_real_)
effects$Importance_pct <- 100*effects$I_k
write_table(effects,file.path(OUT,"process_effects_by_parameter.csv"))
summary <- omnibus <- pairs <- labels <- list()
for(panel_name in unique(interaction(effects$Group,effects$TPL,effects$Condition,drop=TRUE))) {
 d <- effects[interaction(effects$Group,effects$TPL,effects$Condition,drop=TRUE)==panel_name,]
 info <- d[1,c("Group","TPL","Condition")]
 mat <- sapply(PLOT_ORDER,function(p) d$Importance_pct[d$Process==p][order(d$Rank[d$Process==p])])
 if(any(!is.finite(mat))) stop("Undefined importance: all process effects are zero in a parameter set.")
 kw <- if(diff(range(mat))<=ZERO_TOL) list(statistic=0,p.value=1) else kruskal.test(d$Importance_pct,factor(d$Process,levels=PLOT_ORDER))
 pp <- dunn_pairs(d$Importance_pct,d$Process);pp$P_Bonferroni <- p.adjust(pp$P_raw,"bonferroni")
 pp$P_for_letters <- if(kw$p.value>=ALPHA) 1 else pp$P_Bonferroni
 pm <- matrix(1,4,4,dimnames=list(PLOT_ORDER,PLOT_ORDER))
 for(i in seq_len(nrow(pp))) pm[pp$Process1[i],pp$Process2[i]] <- pm[pp$Process2[i],pp$Process1[i]] <- pp$P_for_letters[i]
 ll <- compact_letters(pm,colMeans(mat))
 draws <- replicate(BOOT_REPS,colMeans(mat[sample.int(10,10,replace=TRUE),,drop=FALSE]))
 for(p in PLOT_ORDER) {
   dd <- d[d$Process==p,];q <- quantile(draws[p,],c(.025,.975),names=FALSE)
   summary[[paste(panel_name,p)]] <- cbind(info,data.frame(Process=p,N=10,
      Mean_effect=mean(dd$Effect),SD_effect=sd(dd$Effect),Importance_pct=mean(dd$Importance_pct),
      BOOT95_low=q[1],BOOT95_high=q[2],Positive_fraction=mean(dd$Effect>ZERO_TOL),Letter=ll[p]))
 }
 omnibus[[as.character(panel_name)]] <- cbind(info,Statistic=unname(kw$statistic),P=kw$p.value)
 pairs[[as.character(panel_name)]] <- cbind(info[rep(1,nrow(pp)),],pp)
}
s <- do.call(rbind,summary)
write_table(s,file.path(OUT,"process_importance_summary.csv"))
write_table(do.call(rbind,omnibus),file.path(OUT,"Kruskal_Wallis.csv"))
write_table(do.call(rbind,pairs),file.path(OUT,"Dunn_Bonferroni.csv"))
require_packages("ggplot2")
s$Process <- factor(s$Process,levels=PLOT_ORDER);effects$Process <- factor(effects$Process,levels=PLOT_ORDER)
for(tpl in unique(s$TPL)) {
 d <- s[s$TPL==tpl,];points <- effects[effects$TPL==tpl,]
 p <- ggplot2::ggplot(d,ggplot2::aes(Process,Importance_pct))+ggplot2::geom_col(fill="#8AAFC3")+
   ggplot2::geom_errorbar(ggplot2::aes(ymin=BOOT95_low,ymax=BOOT95_high),width=.2)+
   ggplot2::geom_point(data=points,position=ggplot2::position_jitter(width=.1,seed=731),size=1,alpha=.6)+
   ggplot2::geom_text(ggplot2::aes(y=BOOT95_high+4,label=Letter))+
   ggplot2::facet_grid(Condition~Group)+ggplot2::theme_bw()+
   ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))+
   ggplot2::labs(x=NULL,y="Relative process-removal effect (%)",title=tpl)
 ggplot2::ggsave(file.path(OUT,paste0(tpl,"_process_importance.pdf")),p,width=13,height=6)
}
