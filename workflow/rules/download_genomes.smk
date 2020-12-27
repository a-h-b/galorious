# get taxonomy from GTDB or (if taxonomy isn't checked) from config file

if "taxonomy_check" not in STEPS and config['pangenome_download']['species']:
    localrules: print_species

    rule print_species:
        input:
            OUTPUTDIR
        output:
            "taxonomy/species"
        params:
            spec=config['pangenome_download']['species']
        message: "getting species from config"
        shell:
            """
            echo {params.spec} > {output}
            """

if not config['pangenome_download']['path_to_metadata']:
    rule download_metadata:
        input:
            OUTPUTDIR
        output:
            "tmp/bac120_metadata_r95.tsv"
        resources:
            runtime="4:00:00",
            mem=config['normalMem'] 
        params:
            int="tmp/bac120_metadata_r95.tar.gz",
            url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz"
        threads: 1
        log: "logs/download_metadata.log"
        message: "downloading GTDB metadata"
        shell:
            """
            wget -O {params.int} {params.url}
            tar xvzf {params.int}
            """


if not config['pangenome_download']['path_to_taxonomy']:
    rule download_metadata:
        input:
            OUTPUTDIR
        output:
            "tmp/bac120_taxonomy_r95.tsv"
        resources:
            runtime="4:00:00",
            mem=config['normalMem']
        params:
            int="tmp/bac120_taxonomy_r95.tsv.gz",
            url="https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_taxonomy_r95.tsv.gz"
        threads: 1
        log: "logs/download_taxonomy.log"
        message: "downloading GTDB taxonomy"
        shell:
            """
            wget -O {params.int} {params.url}
            gunzip {params.int}
            """


def input_getDownloadList(step_list,species,metadata,taxonomy):
    if "taxonomy_check" not in STEPS and not species:
        raise Exception("You'd asked to download genomes, but you've not set a species and are not checking the taxonomy.")
    else:
       if metadata:
           inputs=[metadata]
       else:
           inputs=["tmp/bac120_metadata_r95.tsv"]
       if taxonomy:
           inputs.append(taxonomy)
       else:
           inputs.append("tmp/bac120_taxonomy_r95.tsv")
       inputs.append("taxonomy/species")
       return inputs

# if taxonomy_check:
#  use GTDB output to define download:
#   check taxonomy for cluster members of classified taxon
#   if more than max number:
#    use only top max genomes
#   if less than max number:
#    use all
#   if relatives also enabled:
#    check taxonomy for cluster members of related
#    if more than max number:
#     use only reps + top remaining to max relatives' cluster members
#  download list using dataset
# else
#  use GTDB taxonomy to find representatives
#  if perfect match to species, use all matches
#   use representatives + cluster members as far as possible within max
#  if imperfect match, use all imperfect matches (first species + _X, then genus or genus + _X)
#   same thing with representatives
#  ignore relatives 
#  download list using dataset


rule construct_download_list:
    input:
        input_getDownloadList(STEPS,config['pangenome_download']['species'],config['pangenome_download']['path_to_metadata'],
                               config['pangenome_download']['path_to_taxonomy']) 
    output:
        "pangenome/download_list.txt"
    conda:
        ENVDIR + "galorious_roary.yaml"
    threads: 1
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        maxGenome=config['pangenome_download']['max_genomes'],
        useRelatives=config['pangenome_download']['use_gtdb_relatives'],
        steps=STEPS,
        potGTDBresult="taxonomy/GTDB/gtdbtk.bac120.summary.tsv"
    log: "logs/pangenome_download_list.log"
    message: "constructing download list"
    script:
        SCRIPTSDIR + "construct_download_list.R"

checkpoint downloadgenomes:
    input:
        "pangenome/download_list.txt"
    output:
        directory("pangenome/genomes")
    threads: 1
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    conda:
        ENVDIR + "galorious_annotation.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        mkdir -p {output} 
        if grep -q "^GCA" {input}; then
          grep "^GCA" {input} | while read -r acc; do 
            esearch -db assembly -query $acc </dev/null | esummary \
             | xtract -pattern DocumentSummary -element FtpPath_GenBank \
             | while read -r url; do 
                 fname=$(echo $url | grep -o 'GC._.*' | sed 's/$/_genomic.fna.gz/') ;
                 wget -O {output}/$fname "$url/$fname" ;
               done ;
          done
        fi
        if grep -q "^GCF" {input}; then
          grep "^GCF" {input} | while read -r acc; do 
            esearch -db assembly -query $acc </dev/null | esummary \
             | xtract -pattern DocumentSummary -element FtpPath_RefSeq \
             | while read -r url; do 
                 fname=$(echo $url | grep -o 'GC._.*' | sed 's/$/_genomic.fna.gz/') ;
                 wget -O {output}/$fname "$url/$fname" ;
               done ;
          done
        fi
        """



localrules: cp_genomes

rule cp_genomes:
    input:
        "pangenome/genomes"
    output:
        touch("status/download_pangenome.done") 

def gather_genomes(wildcards):
    checkpoint_output=checkpoints.downloadgenomes.get().output[0]
    return expand("pangenome/annotation/{i}.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)


