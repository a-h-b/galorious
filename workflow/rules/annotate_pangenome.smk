rule prokkaC_references:
    input:
        "pangenome/genomes/{genome}.fna.gz",
        "status/download_pangenome.done" 
    output:
        "pangenome/annotation/{genome}.gff",
        "pangenome/annotation/{genome}.faa",
        "pangenome/annotation/{genome}.fna",
        "pangenome/annotation/{genome}.ffn",
        "pangenome/annotation/{genome}.fsa"
    threads: 1 
#       getThreads(8)
    log: "logs/pangenomes_annotate.{genome}.log"
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        out="pangenome/annotation"
    conda:
        ENVDIR + "galorious_annotation.yaml"
    message: "prokkaC_references: annotating {wildcards.genome}."
    shell:
        """
        IN={input[0]}
        IN_UC=${{IN%.*}}
        zcat {input[0]} > $IN_UC
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then
          {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --setupdb
        fi
        {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --force --outdir {params.out} --prefix {wildcards.genome} --noanno --cpus {threads} $IN_UC >> {log} 2>&1
        rm $IN_UC
        """


if not 'download_pangenome' in STEPS and config['inputs']['Genomes2Compare']:
    localrules: cp_genome
    checkpoint cp_genomes:
        input:
            config['inputs']['Genomes2Compare']
        output:
            directory("pangenome/genomes")
        shell:
            """
            mkdir -p {output}
            cp {input}/*fna.gz {output}
            """

    def gather_genomes(wildcards):
        checkpoint_output=checkpoints.cp_genomes.get().output[0]
        return expand("pangenome/annotation/{i}.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i) 

localrules: ctrl_annotate_pangenome
if not config['inputs']['Gffs2Compare']:
    rule ctrl_annotate_pangenome:
        input:
            gather_genomes
        output:
            touch("status/annotate_pangenome.done")
else:
    rule ctrl_annotate_pangenome:
        input:
            config['inputs']['Gffs2Compare'],
            gather_genomes
        output:
            "status/annotate_pangenome.done"
        params:
            annodir="pangenome/annotation"
        shell:
            """
            cp {input[0]}/*gff {params.annodir}/ && touch {output}
            """
