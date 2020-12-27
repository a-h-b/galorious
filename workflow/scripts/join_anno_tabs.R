log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

pafile <- snakemake@input[[1]]
hmmfiles <- unlist(snakemake@input[-1])
pa_outfile <- snakemake@output[[1]]

pa <- read.delim(pafile,stringsAsFactors=F,sep=",")
first <- T
for(f in hmmfiles){
  if(first){
    annos <- read.delim(f,stringsAsFactors=F)
    canno <- colnames(annos)[2]
    colnames(annos)[2] <- "annotation"
    annos$type <- canno
    first <- F
  }else{
    tmpa <- read.delim(f,stringsAsFactors=F)
    canno <- colnames(tmpa)[2]
    colnames(tmpa)[2] <- "annotation"
    tmpa$type <- canno
    annos <- rbind(annos,tmpa)
  }
}
anno <- as.data.frame(tapply(annos$annotation,list(annos$group,annos$type),function(x) x),
         stringsAsFactors=F)
pa <- merge(pa,anno,by.x=1,by.y=0,all.x=T)

saveRDS(pa,pa_outfile)


