rule prokkaC_references:
    input:
        "pangenome/genomes/{genome}.fna.gz",
        "status/download_pangenome.done" 
    output:
        "pangenome/annotation/{genome}.prokka.gff",
        "pangenome/annotation/{genome}.prokka.faa",
        "pangenome/annotation/{genome}.prokka.fna",
        "pangenome/annotation/{genome}.prokka.ffn",
        "pangenome/annotation/{genome}.prokka.fsa"
    threads: getThreads(8)
    log: "logs/pangenomes_annotate.{genome}.log"
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        out="pangenome/annotation",
        pref="{wildcards.genome}.prokka"
    conda:
        ENVDIR + "galorious_annotation.yaml"
    message: "prokkaC_references: annotating {wildcards.genome}."
    shell:
        """
        IN={input[0]}
        IN_UC=${{IN%.*}}
        zcat {input} > $IN_UC
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then
          {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --setupdb
        fi
        {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --force --outdir {params.out} --prefix {params.pref} --noanno --cpus {threads} $IN_UC >> {log} 2>&1
        rm $IN_UC
        """


def gather_genomes(wildcards):
    checkpoint_output=checkpoints.cp_genomes.get().output[0]
    return expand("pangenome/annotation/{i}.prokka.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.fna.gz")).i)

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
