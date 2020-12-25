def getFocalGff(step_list,Gff):
    if "annotation" in step_list:
        return ["pangenome/annotation","annotation/prokka.gff"]
    elif Gff:
        return ["pangenome/annotation",Gff]
    else:
        return "pangenome/annotation"

def gather_annos(step_list,custom_genomes):
    if 'annotate_pangenome' in step_list:
        checkpoint_output=checkpoints.downloadgenomes.get(**wildcards).output[0]
        all_annos = expand("pangenome/annotation/{i}.prokka.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)
    elif custom_genome:
        checkpoint_output=checkpoints.cpgenomes.get(**wildcards).output[0]
        all_annos = expand("pangenome/annotation/{i}.prokka.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)


if not 'annotate_pangenome' in step_list and config['inputs']['Gffs2Compare']:
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
        directory("pangenome/roary"),
        "status/pangenomics.done"
    threads:
        getThreads(8)
    resources:
        runtime="48:00:00",
        mem=config['normalMem']
    log: "logs/roary.log"
    conda:
        ENVDIR + "galorious_roary.yaml"
    shell:
        """
        roary -p {threads} -f {output[0]} -e --mafft {input} &>> {log} 
        touch {output[1]}
        """
