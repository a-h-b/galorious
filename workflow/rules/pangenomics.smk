def getFocalGff(step_list,Gff):
    all = []
    if "annotation" in step_list:
        all.append("annotation/prokka.gff")
    if "annotate_pangenome" in step_list:
        all.append(gather_genomes)
    elif GFF:
        all.append("pangenome/annotation")
    return all

#def gather_annos(step_list,custom_genomes):
#    if 'annotate_pangenome' in step_list:
#        checkpoint_output=checkpoints.downloadgenomes.get(**wildcards).output[0]
#        all_annos = expand("pangenome/annotation/{i}.prokka.gff",
#                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)
#    elif custom_genome:
#        checkpoint_output=checkpoints.cp_genomes.get(**wildcards).output[0]
#        all_annos = expand("pangenome/annotation/{i}.prokka.gff",
#                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)


if not 'annotate_pangenome' in STEPS and config['inputs']['Gffs2Compare']:
    localrules: cp_gffs
    checkpoint cp_gffs:
        input:
            config['inputs']['Gffs2Compare']
        output:
            directory("pangenome/annotation")
        shell:
            """
            mkdir -p {output}
            cp {input}/*.gff {output}
            """



rule roary:
    input:
        getFocalGff(STEPS,config['inputs']['Gff'])
    output:
        "status/pangenomics.done",
        "pangenome/roary/pan_genome_reference.fa",
        "pangenome/roary/clustered_proteins",
        "pangenome/roary/gene_presence_absence.csv",
        "pangenome/roary/core_accessory_graph.dot"
    threads: 1
#        getThreads(8)
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    params:
        outdir="pangenome/roary"
    log: "logs/roary.log"
    conda:
        ENVDIR + "galorious_roary.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        roary -p {threads} -f {params.outdir} -e --mafft {input} &>> {log} 
        touch {output[0]}
        """

rule convert_roary_nets:
    input:
        "pangenome/roary/core_accessory_graph.dot"
    output:
        "pangenome/roary/core_accessory_graph.RDS"
    threads: 1
    resources:
        runtime="1:00:00",
        mem=config['normalMem']
    log: "logs/roary_network.log"
    conda:
        ENVDIR + "galorious_roary.yaml"
    script:
        SCRIPTSDIR + "dot2edgelist.R"


rule annotate_roary_focal:
    input:
        "pangenome/roary/gene_presence_absence.csv",
        "pangenome/roary/clustered_proteins"
    output:
        "pangenome/roary/focal_gene_presence_absence.Rtab"
        threads: 1
    resources:
        runtime="2:00:00",
        mem=config['normalMem']
    log: "logs/roary_focal_presence.log"
    conda:
        ENVDIR + "galorious_roary.yaml"
    script:
        SCRIPTSDIR + "focal_gene_presence.R"
