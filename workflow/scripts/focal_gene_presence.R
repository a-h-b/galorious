log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

dot2edge <- function(dotstring){
  edges <- data.frame("string"=dotstring[(grep("[*]",dotstring)[2]+1):
                                         (grep("[}]",dotstring[,1])-1)],
                      stringsAsFactors=F)
  edges$source <- gsub("^ *","",gsub(" ->.+","",edges$string))
  edges$sink <- gsub(" .+","",gsub(".+-> ","",edges$string))
  edges$weight <- as.numeric(gsub(".+weight = ","",gsub("];$","",edges$string)))
  edges[,c("source","sink","weight")] # to be a full edgelist, the edgelist has to be duplicated with source/sink reversed, but it's not done here 
}

pafile <- snakemake@input[[1]]
clusterfile <- snakemake@input[[2]]
pa_outfile <- snakemake@output[[1]]

pa <- read.delim(pafile,stringsAsFactors=F)
saveRDS(dot2edge(dotread[,1]),edgefile)


