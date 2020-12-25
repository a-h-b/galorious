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
           inputs=metadata
       else:
           inputs="tmp/bac120_metadata_r95.tsv"
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

rule downloadgenomes:
    input:
        "pangenome/download_list.txt"
    output:
        directory("pangenome/genomes")
    threads: 1
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    conda:
        ENVDIR + "galorious_roary.yaml"
    shell:
        """
        if [ -f $CONDA_PREFIX/bin/datasets ]; then
           curl -o $CONDA_PREFIX/bin/datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
        fi
        chmod +x $CONDA_PREFIX/bin/datasets
        mkdir -p {output} && cd {output}
        datasets download genome accession --inputfile ../../{input[0]} --exclude-gff3 --exclude-protein --exclude-rna
        """

localrules: cpgenomes

if config['inputs']['Genomes2Compare']:
    checkpoint cpgenomes:
        input:
            "pangenome/genomes",
            config['inputs']['Genomes2Compare']
        output:
            "status/download_pangenome.done" 
        shell:
            """
            cp {input[1]}/*fna.gz {input[0]} && touch {output}
            """
else:
    checkpoint cpgenomes:
        input:
            "pangenome/genomes"
        output:
            touch("status/download_pangenome.done") 

