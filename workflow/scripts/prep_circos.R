contigs <- read.delim("assembly/assembly.length.txt",header=F,stringsAsFactors=F)
colnames(contigs) <- c("contig","length")
genes <- read.delim("annotation/annotation_CDS_RNA_hmms.gff",header=F,stringsAsFactors=F)
colnames(genes) <- c("contig","source","feature","start","end","score","strand","frame","attribute")
skew <- read.delim("assembly/assembly.gc_skew.txt",stringsAsFactors=F)
covSelf <- read.delim("assembly/assembly.contig_depth_perBase.txt",stringsAsFactors=F,header=F)
colnames(covSelf) <- c("contig","pos","depth")

dir.create("visualization/data", showWarnings = FALSE,recursive = T) 

kt <- data.frame("chr"=rep("chr",nrow(contigs)),
                   "ps"=rep("-",nrow(contigs)),
                   "ID"=contigs$contig,"LABEL"=contigs$contig,
                   "START"=rep(0,nrow(contigs)),"END"=contigs$length,
                   "COLOR"=rep("black",nrow(contigs)),stringsAsFactors=F)
kt <- kt[order(kt$END,decreasing=T),]
write.table(kt,"visualization/data/contigs.ideogram",row.names=F,col.names=F,quote=F)
  

kt <- data.frame("chr"=genes$contig[genes$strand=="+"&genes$feature=="CDS"],
                   "START"=genes$start[genes$strand=="+"&genes$feature=="CDS"],
                   "END"=genes$end[genes$strand=="+"&genes$feature=="CDS"],
                   "COLOR"=paste("color=",
                                 ifelse(grepl("partial=00",genes$attribute[genes$strand=="+"&genes$feature=="CDS"]),
                                                 "(51,51,51)","(179,179,179)"),sep=""),
                   stringsAsFactors=F)
write.table(kt,"visualization/data/features.forward.txt",row.names=F,col.names=F,quote=F)
kt <- data.frame("chr"=genes$contig[genes$strand=="-"&genes$feature=="CDS"],
                   "START"=genes$start[genes$strand=="-"&genes$feature=="CDS"],
                   "END"=genes$end[genes$strand=="-"&genes$feature=="CDS"],
                   "COLOR"=paste("color=",ifelse(grepl("partial=00",genes$attribute[genes$strand=="-"&genes$feature=="CDS"]),
                                                      "(51,51,51)","(179,179,179)"),sep=""),
                   stringsAsFactors=F)
write.table(kt,"visualization/data/features.reverse.txt",row.names=F,col.names=F,quote=F)

if(length(genes$contig[genes$feature=="rRNA"])>0){
  kt <- data.frame("chr"=genes$contig[genes$feature=="rRNA"],
                   "START"=genes$start[genes$feature=="rRNA"],
                   "END"=genes$end[genes$feature=="rRNA"],
                   "COLOR"=paste("color=",ifelse(genes$strand[genes$feature=="rRNA"]=="+",
                                                 "vvdgreen","vvdblue"),sep=""),
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/rRNAs.txt",row.names=F,col.names=F,quote=F)

if(length(genes$contig[genes$strand=="+"&genes$feature=="region"])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="+"&genes$feature=="region"],
                   "START"=genes$start[genes$strand=="+"&genes$feature=="region"],
                   "END"=genes$end[genes$strand=="+"&genes$feature=="region"],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/regions.forward.txt",row.names=F,col.names=F,quote=F)
if(length(genes$contig[genes$strand=="-"&genes$feature=="region"])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="-"&genes$feature=="region"],
                   "START"=genes$start[genes$strand=="-"&genes$feature=="region"],
                   "END"=genes$end[genes$strand=="-"&genes$feature=="region"],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/regions.reverse.txt",row.names=F,col.names=F,quote=F)

if(length(genes$contig[genes$strand=="+"&
                       (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                        grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                        grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="+"&
                                        (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                           grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                           grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   "START"=genes$start[genes$strand=="+"&
                                         (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                            grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                            grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   "END"=genes$end[genes$strand=="+"&
                                     (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                        grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                        grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/goi.forward.txt",row.names=F,col.names=F,quote=F)
if(length(genes$contig[genes$strand=="-"&genes$feature=="region"])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="-"&
                                        (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                           grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                           grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   "START"=genes$start[genes$strand=="-"&
                                         (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                            grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                            grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   "END"=genes$end[genes$strand=="-"&
                                     (grepl("KEGG=K01505",genes$attribute)|grepl("Pfam_A=BCCT",genes$attribute)|
                                        grepl("Pfam_A=OpuAC",genes$attribute)|grepl("dbCAN=GH43",genes$attribute)|
                                        grepl("dbCAN=GH29.hmm",genes$attribute)|grepl("dbCAN=GH10.hmm",genes$attribute))],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/goi.reverse.txt",row.names=F,col.names=F,quote=F)

if(length(genes$contig[genes$strand=="+"&
                       (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                        grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                        grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                        grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                        grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="+"&
                                        (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                           grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                           grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                           grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                           grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   "START"=genes$start[genes$strand=="+"&
                                         (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                            grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                            grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                            grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                            grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   "END"=genes$end[genes$strand=="+"&
                                     (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                        grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                        grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                        grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                        grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/motu.forward.txt",row.names=F,col.names=F,quote=F)
if(length(genes$contig[genes$strand=="-"&genes$feature=="region"])>0){
  kt <- data.frame("chr"=genes$contig[genes$strand=="-"&
                                        (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                           grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                           grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                           grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                           grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   "START"=genes$start[genes$strand=="-"&
                                         (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                            grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                            grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                            grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                            grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   "END"=genes$end[genes$strand=="-"&
                                     (grepl("essential=TIGR00092",genes$attribute)|grepl("essential=tRNA-synt_1d",genes$attribute)|
                                        grepl("essential=TIGR00435",genes$attribute)|grepl("essential=TIGR00959",genes$attribute)|
                                        grepl("essential=TIGR00422",genes$attribute)|grepl("essential=TIGR00468",genes$attribute)|
                                        grepl("essential=TIGR00414",genes$attribute)|grepl("essential=TIGR00396",genes$attribute)|
                                        grepl("KEGG=K01409",genes$attribute)|grepl("essential=TIGR00064",genes$attribute))],
                   stringsAsFactors=F)
} else {
  kt <- ""
}
write.table(kt,"visualization/data/motu.reverse.txt",row.names=F,col.names=F,quote=F)



kt <- data.frame("chr"=skew$Sequence,
                 "START"=skew$Index,
                 "END"=skew$Index+1999,
                 "VALUE"=skew$GC.Skew..2kb.,
                 stringsAsFactors=F)
for(i in unique(kt$chr)){
  kt$END[kt$chr==i][which.max(kt$END[kt$chr==i])] <- contigs$length[contigs$contig==i]
}
options("scipen"=100, "digits"=4)
write.table(kt,"visualization/data/gc_skew.txt",row.names=F,col.names=F,quote=F)

kt$VALUE <- sapply(1:nrow(kt),
                   function(x) mean(covSelf$depth[covSelf$contig==kt$chr[x]&
                                                    covSelf$pos>=kt$START[x]&
                                                    covSelf$pos<=kt$END[x]]))
write.table(kt,"visualization/data/cov_self.txt",row.names=F,col.names=F,quote=F)

for(g in list.files(path="additional_mapping",pattern="contig_depth_perBase.txt")){
  sample <- gsub(".contig_depth_perBase.txt","",g)
  covA <- read.delim(paste0("additional_mapping/",g),stringsAsFactors=F,header=F)
  colnames(covA) <- c("contig","pos","depth")
  kt$VALUE <- sapply(1:nrow(kt),
                     function(x) mean(covA$depth[covA$contig==kt$chr[x]&
                                                      covA$pos>=kt$START[x]&
                                                      covA$pos<=kt$END[x]]))
  write.table(kt,paste0("visualization/data/cov_",sample,".txt"),
              row.names=F,col.names=F,quote=F)
}

if(file.exists("pangenome/roary/gene_presence_absence.anno.RDS")){
  pan <-readRDS("pangenome/roary/gene_presence_absence.anno.RDS")
  kt <- data.frame("chr"=genes$contig[genes$feature=="CDS"],
                   "START"=genes$start[genes$feature=="CDS"],
                   "END"=genes$end[genes$feature=="CDS"],
                   "VALUE"=abs(sapply(gsub("ID=","",gsub(";.+","",
                                                         genes$attribute[genes$feature=="CDS"])),
                                      function(x){
                      if(any(pan$prokka==x)) pan$No..isolates[pan$prokka==x] else max(pan$No..isolates)
                   })-max(pan$No..isolates)),
                   stringsAsFactors=F)
}else {
  kt <- ""
}
write.table(kt,"visualization/data/pangenome.txt",
            row.names=F,col.names=F,quote=F)



