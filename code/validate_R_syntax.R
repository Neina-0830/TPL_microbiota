# Parse every R script without running the analyses.
root <- normalizePath(Sys.getenv("TPL_PROJECT_ROOT",unset=getwd()),mustWork=TRUE)
files <- list.files(file.path(root,"code"),pattern="\\.R$",recursive=TRUE,full.names=TRUE)
for(f in files) parse(file=f)
cat("PASS:",length(files),"R files parsed.\n")
