# Install the dependencies used by the analysis modules.
pkgs <- c("ape","bigmemory","dplyr","emmeans","ggplot2","iCAMP",
          "lme4","lmerTest","metafor","performance","pbkrtest",
          "piecewiseSEM","tibble","tidyr","vegan","permute")
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) install.packages(missing,repos="https://cloud.r-project.org")
if(any(!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE))) stop("Some dependencies could not be installed.")
dir.create("results",showWarnings=FALSE)
capture.output(sessionInfo(),file="results/R_sessionInfo.txt")
