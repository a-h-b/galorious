def getFocalGff(step_list,Gff):
    all = []
    if "annotation" in step_list:
        all.append("annotation/prokka.gff")
    if "annotate_pangenome" in step_list:
        all.append(gather_genomes)
    elif GFF:
        all.append("pangenome/annotation")
    return all


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

rule fastANI:
    input:
        "assembly/unicycler/assembly.fasta",
        "pangenome/genomes"
    output:
        "pangenome/fastANI/all.tsv"
    params:
        outdir="pangenome/fastANI"
    threads: 1
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    conda:
        ENVDIR + "galorious_gtdbtk.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        for genome in `ls {input[1]}/*fna.gz`
        do 
          OUT=$(basename $genome .gz)
          fastANI -r {input[0]} -q $genome -o {params.outdir}/ani.$OUT.txt  --visualize
        done
        cat {params.outdir}/ani.$OUT.txt >> {output} 
        """

rule roary:
    input:
        getFocalGff(STEPS,config['inputs']['Gff'])
    output:
        "pangenome/roary/pan_genome_reference.fa",
        "pangenome/roary/clustered_proteins",
        "pangenome/roary/gene_presence_absence.csv",
        "pangenome/roary/core_accessory_graph.dot"
    threads: 4
#        getThreads(8)
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        outdir="pangenome/roary",
        tmp="pangenome/roary_"
    log: "logs/roary.log"
    conda:
        ENVDIR + "galorious_roary.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        roary -z -p {threads} -f {params.outdir} -e --mafft {input} &>> {log}
        mv {params.tmp}*/* {params.outdir}/
        rm -r {params.tmp}*
        tar cvzf {params.outdir}/intermediary.tar.gz {params.outdir}/_*
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


rule translate_pangenome:
    input:
        "pangenome/roary/pan_genome_reference.fa"
    output:
        "pangenome/roary/pan_genome_reference.faa"
    threads: 1
    resources:
        runtime="1:00:00",
        mem=config['normalMem']
    log: "logs/roary_network.log"
    conda:
        ENVDIR + "galorious_utils.yaml"
    shell:
        """
        transeq -table 11 {input} {output} &>> {log}
        sed -i 's/>.*group/>group/' {output}
        """

rule pangenomics_hmmer:
    input:
        "pangenome/roary/pan_genome_reference.faa"
    output:
        "pangenome/annotation/pan_genome_reference.faa.{db}.hmmscan"
    params:
        lambda wildcards: config["hmm_settings"][wildcards.db]["cutoff"],
        dbs = DBPATH + "/hmm/{db}"
    resources:
        runtime = "6:00:00",
        mem = config['normalMem']
    conda: ENVDIR + "galorious_annotation.yaml"
    threads: getThreads(4)
    log: "logs/analysis_hmmer.{db}.log"
    message: "hmmer: Running HMMER for {wildcards.db}."
    shell:
        """
        sed 's/>.* group/>group/' {input} {input}
        hmmsearch --cpu {threads} {params[0]} --noali --notextw \
          --tblout {output} {params.dbs}/*.hmm {input} >/dev/null 2> {log}
        """

rule pangenomics_hmm2tab:
    input:
        "pangenome/annotation/pan_genome_reference.faa.{db}.hmmscan"
    output:
        "pangenome/annotation/pangenome.anno.{db}.tsv"
    resources:
        runtime = "4:00:00",
        mem = config['normalMem']
    threads: 1
    params:
        "{db}",
        lambda wildcards: config["hmm_settings"][wildcards.db]["trim"],
        "pangenome/roary/gene_presence_absence.csv"
    log: "logs/pangenomics_hmm2tab.{db}.log"
    message: "hmm2tab: Best hmmer result for {wildcards.db}."
    conda:
        ENVDIR + "galorious_annotation.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SCRIPTSDIR}/hmmscan_addBest2tab.pl -i {input} -a {params[2]} \
         -n {params[0]} {params[1]} -o {output} -g $(wc -l {params[2]}) > {log} 2>&1
        """

rule merge_hmmtabs:
    input:
        "pangenome/roary/gene_presence_absence.csv",
        expand("pangenome/annotation/pangenome.anno.{db}.tsv",db=config["hmm_DBs"].split())
    output:
        "pangenome/roary/gene_presence_absence.anno.RDS"
    resources:
        runtime = "4:00:00",
        mem = config['normalMem']
    threads: 1
    log: "logs/pangenomics_mergehmm.log"
    message: "hmm2tab: Best hmmer result per database."
    conda: 
        ENVDIR + "galorious_roary.yaml"
    script:
        SCRIPTSDIR + "join_anno_tabs.R"


localrules: ctrl_pangenomics

rule ctrl_pangenomics:
    input:
        "pangenome/fastANI/all.tsv",
        "pangenome/roary/gene_presence_absence.anno.RDS",
        "pangenome/roary/core_accessory_graph.RDS"
    output:
        touch("status/pangenomics.done")


