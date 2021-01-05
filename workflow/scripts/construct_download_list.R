log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

condap <- Sys.getenv("CONDA_PREFIX")
.libPaths(paste0(condap,"/lib/R/library"))

limit_genome_list <- function(maxNo, meta=metadata){
  meta <- meta[,c("gtdb_type_species_of_genus","gtdb_type_designation","ncbi_type_material_designation",
                  "ncbi_genome_representation","checkm_completeness","checkm_contamination","ncbi_genome_category",
                  "contig_count","ambiguous_bases","ssu_length","accession","gtdb_genome_representative")]
 meta$score <- 100000*as.numeric(meta$gtdb_type_species_of_genus=="t"&!is.na(meta$gtdb_type_species_of_genus)) + 
               10000*as.numeric(meta$ncbi_type_material_designation!="none"&!is.na(meta$ncbi_type_material_designation)) +
               1000*as.numeric(meta$gtdb_type_designation!="not type material"&!is.na(meta$gtdb_type_designation)) +
               100*as.numeric(meta$ncbi_genome_representation=="full"&!is.na(meta$ncbi_genome_representation)) + 
               ifelse(is.na(as.numeric(meta$checkm_completeness)),0,as.numeric(meta$checkm_completeness)) +
               (ifelse(is.na(as.numeric(meta$checkm_contamination)),0,as.numeric(meta$checkm_contamination)) * -5) +
               (-100 * as.numeric(meta$ncbi_genome_category != "none"&!is.na(meta$ncbi_genome_category))) +
               (-5 * (ifelse(is.na(as.numeric(meta$ambiguous_bases)),0,as.numeric(meta$ambiguous_bases)) / 10000)) +
               10 * as.numeric(as.numeric(meta$ssu_length)>1300&!is.na(meta$ssu_length))
 meta$accession[order(meta$score,decreasing=T)[1:maxNo]]
}

metadata_file <- snakemake@input[[1]]
taxonomy_file <- snakemake@input[[2]]
species_file  <- snakemake@input[[3]]
gtdb_out_file <- snakemake@params[["potGTDBresult"]]

output_file <- snakemake@output[[1]]

maxGenomes <- as.numeric(snakemake@params[["maxGenome"]])
useRelatives <- snakemake@params[["useRelatives"]]
steps <- unlist(snakemake@params[["steps"]])


taxString <- read.delim(species_file,header=F,stringsAsFactors=F)[1,1]
taxon <- as.character(read.delim(species_file,stringsAsFactors=F,header=F,sep=";")[1,])
names(taxon) <- c("d","p","c","o","f","g","s")
taxon <- gsub("^.__","",taxon)

taxonomy <- read.delim(taxonomy_file, header=F,stringsAsFactors=F)
colnames(taxonomy) <- c("ref_genome","taxString")

metadata <- read.delim(metadata_file,stringsAsFactors=F)


if(any(grepl("taxonomy_check",steps)) | taxString %in% taxonomy$taxString){
  focal_ref <- taxonomy$ref_genome[taxonomy$taxString==taxString]
}
if(!"focal_ref" %in% ls() | length(focal_ref)==0){
  if(taxon["s"]!="" & any(grepl(paste0("s__",taxon["s"]), gsub("_. ","",gsub("_.$","",taxonomy$taxString))))){
    focal_ref <- taxonomy$ref_genome[which(grepl(paste0("s__",taxon["s"]), gsub("_. ","",gsub("_.$","",taxonomy$taxString))))]
  }else{
    if(taxon["g"]!="" & any(grepl(paste0("g__",taxon["g"]), gsub("_. ","",gsub("_.$","",taxonomy$taxString))))){
     focal_ref <- taxonomy$ref_genome[which(grepl(paste0("g__",taxon["g"]), gsub("_. ","",gsub("_.$","",taxonomy$taxString))))]
    }
  } 
}
focal_ref <- unique(focal_ref)
if(length(focal_ref)==0){
  print("no reference genome found")
  all <- focal_ref
}else{
  if(length(focal_ref)>maxGenomes){
    focal_ref <- limit_genome_list(maxGenomes, metadata[metadata$accession %in% focal_ref,])
    all <- focal_ref
  }else{
    focal_cluster <- unique(metadata$accession[metadata$gtdb_genome_representative %in% focal_ref &! metadata$accession %in% focal_ref])
    if(length(focal_cluster)+ length(focal_ref) > maxGenomes ){
      focal_cluster <- limit_genome_list(maxGenomes - length(focal_ref),metadata[metadata$accession %in% focal_cluster,])
      all <- c(focal_ref,focal_cluster)
    }else{
      if(useRelatives & any(grepl("taxonomy_check",steps))){
        allGTDB <- read.delim(gtdb_out_file,stringsAsFactors=F)
        relatives <- gsub(" ","",gsub(",.+","",unlist(strsplit(unique(allGTDB$other_related_references.genome_id.species_name.radius.ANI.AF.,split=";")))))
        if(length(focal_cluster) + length(relatives) + length(focal_ref) > maxGenomes){
          rel_ANI <- sapply(unlist(strsplit(unique(allGTDB$other_related_references.genome_id.species_name.radius.ANI.AF.,split=";"))),
                             function(x) as.numeric(gsub(" ","",unlist(strsplit(x,split=","))[4])))
          relatives <- relatives[order(rel_ANI,decreasing=T)[1:(maxGenomes-length(focal_ref)-length(focal_cluster))]]
          all <- c(focal_ref,focal_cluster,relatives)
        }else{
          rel_cluster <- unique(metadata$accession[metadata$gtdb_genome_representative %in% relatives &! metadata$accession %in% relatives])
          if(length(rel_cluster) + length(focal_cluster) + length(relatives) + length(focal_ref) > maxGenomes) rel_cluster <- limit_genome_list(maxGenomes - length(focal_ref) -length(focal_cluster)-length(relatives),metadata[metadata$accession %in% rel_cluster,])
          all <- c(focal_ref,focal_cluster,relatives,rel_cluster)
        }       
      }else{
        all <- c(focal_ref,focal_cluster)
      } 
    }
  } 
}


if(length(all)>1){
  all <- gsub("^.._","",all)
  if(!dir.exists(dirname(output_file))) dir.create(dirname(output_file),recursive =T)
  write.table(data.frame("ass"=unique(all),stringsAsFactors=F),
                                    output_file,
                                    sep="\t",quote=F,row.names=F,col.names=F) 
}

