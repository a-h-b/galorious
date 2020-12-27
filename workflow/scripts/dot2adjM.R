log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

dot2adjM <- function(dotstring){
  edges <- data.frame("string"=dotstring[(grep("[*]",dotstring)[2]+1):
                                         (grep("[}]",dotstring[,1])-1)],
                      stringsAsFactors=F)
  edges$source <- gsub("^ *","",gsub(" ->.+","",edges$string))
  edges$sink <- gsub(" .+","",gsub(".+-> ","",edges$string))
  edges$weight <- as.numeric(gsub(".+weight = ","",gsub("];$","",edges$string)))
  nodes <- dotstring[(grep("[*]",dotstring)[1]+1):(grep("[*]",dotstring)[2]-1)]
  adj <- matrix(0,nrow=length(nodes),ncol=length(nodes),
                dimnames=list(nodes,nodes))
  for(i in 1:nrow(edges)){
    adj[rownames(adj)==edges$source[i],colnames(adj)==edges$sink[i]] <- 1
    adj[rownames(adj)==edges$sink[i],colnames(adj)==edges$source[i]] <- 1   # despite the error, the network isn't directed.  
  }
  adj
}

dotfile <- snakemake@input[[1]]
adjMfile <- snakemake@output[[1]]

dotread <- read.delim(dotfile,stringsAsFactors=F)
saveRDS(dot2adjM(dotread[,1]),adjMfile)


