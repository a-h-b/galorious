rule prokkaC:
    input:
        "pangenome/assemblies/{refseq}_genomic.fna.gz 
    output:
        "pangenome/annotation/{refseq}_genomic.prokka.gff",
        "pangenome/annotation/{refseq}_genomic.prokka.faa",
        "pangenome/annotation/{refseq}_genomic.prokka.fna",
        "pangenome/annotation/{refseq}_genomic.prokka.ffn",
        "pangenome/annotation/{refseq}_genomic.prokka.fsa"
    threads: getThreads(8)
    log: "logs/pangenomes_annotate.{refseq}.log"
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        out="pangenome/annotation",
        pref="{wildcards.refseq}_genomic.prokka"
    conda:
        ENVDIR + "galorious_annotation.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        if [ ! -f $CONDA_PREFIX/db/hmm/HAMAP.hmm.h3m ]; then
          {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --setupdb
        fi
        {BINDIR}/prokkaC --dbdir $CONDA_PREFIX/db --force --outdir {params.out} --prefix {params.pref} --noanno --cpus {threads} {input[0]} >> {log} 2>&1
        """

# get taxonomy from GTDB or (if taxonomy isn't checked) from config file

checkpoint dry_ncbi_all:
    input:
        "taxonomy/GTDB/gtdbtk.bac120.summary.tsv"
    output:
        "pangenome/ncbi_dry_results.txt"
    conda:
        "ENVDIR + "galorious_roary.yaml"
    threads: 1
    message: "checking ncbi for genomes"
    shell:
        """
        TAX=`cut -f {input}   | tail -n 1 | sed 's#.*g__##' | cut -f 1 -d ";"`
        ncbi-genome-download --formats fasta --genera $TAX --dry-run bacteria | wc -l >> {output}
        """

checkpoint ncbi_all:
    input:
        "taxonomy/GTDB/gtdbtk.bac120.summary.tsv"
    output:
        genomes=directory("pangenome/ncbi_all/refseq/bacteria")
    params:
        intdir="pangenome/ncbi_all"
    conda:
        "ENVDIR + "galorious_roary.yaml"
    threads: 1
    message: "getting genomes from ncbi"
    shell:
        """
        TAX=`cut -f {input}   | tail -n 1 | sed 's#.*g__##' | cut -f 1 -d ";"`
        mkdir -p {params.intdir} && cd {params.intdir}
        ncbi-genome-download --formats fasta --genera $TAX bacteria 
        """

checkpoint ncbi_representatives:
    input:
        "taxonomy/GTDB/gtdbtk.bac120.summary.tsv"
    output:
        genomes=directory("pangenome/ncbi_representatives/refseq/bacteria")
    params:
        intdir="pangenome/ncbi_representatives"
    conda:
        "ENVDIR + "galorious_roary.yaml"
    threads: 1
    message: "getting representative genomes from ncbi"
    shell:
        """
        TAX=`cut -f {input}   | tail -n 1 | sed 's#.*g__##' | cut -f 1 -d ";"`
        mkdir -p {output} && cd {output}
        ncbi-genome-download --formats fasta --genera $TAX  --refseq-categories reference bacteria      
        """

# input function for the rule aggregate
def aggregate_dryrun():
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with checkpoints.dry_ncbi_all.get().output[0].open() as f:
        if int(f.read().strip()) > 200:
            return "pangenome/ncbi_representatives"
        else:
            return "pangenome/ncbi_all"


def aggregate_downloaded(wildcards):
    checkpoint_output_1 = checkpoints.dry_ncbi_all.get().output[0]
    with checkpoints.dry_ncbi_all.get().output[0].open() as f:
        if int(f.read().strip()) > 200:
            checkpoint_output_2 = checkpoints.ncbi_representatives.get().output[0]
            return expand("pangenome/ncbi_representatives/refseq/bacteria/{i}.fna.gz",
                          i=glob_wildcards(os.path.join(checkpoint_output_2, "{i}.fna.gz")).i)
        else:
            checkpoint_output_2 = checkpoints.ncbi_all.get().output[0]
            return expand("pangenome/ncbi_all/refseq/bacteria/{i}.fna.gz",
                          i=glob_wildcards(os.path.join(checkpoint_output_2, "{i}.fna.gz")).i)

def aggregate_forProkka(wildcards):
    checkpoint_output_1 = checkpoints.dry_ncbi_all.get().output[0]
    with checkpoints.dry_ncbi_all.get().output[0].open() as f:
        if int(f.read().strip()) > 200:
            checkpoint_output_2 = checkpoints.ncbi_representatives.get().output[0]
            return expand("epresentatives/refseq/bacteria/{i}.fna.gz",
                          i=glob_wildcards(os.path.join(checkpoint_output_2, "{i}.fna.gz")).i)
        else:
            checkpoint_output_2 = checkpoints.ncbi_all.get().output[0]
            return expand("pangenome/ncbi_all/refseq/bacteria/{i}.fna.gz",
                          i=glob_wildcards(os.path.join(checkpoint_output_2, "{i}.fna.gz")).i)



localrules: aggregate_ncbi,aggregate_refseq_genomes

rule aggregate_ncbi:
    input:
        aggregate_dryrun
    output:
        "pangenome/assemblies"
    shell:
        "mkdir -p {output}"

rule aggregate_refseq_genomes:
    input:
        aggregate_dryrun
    output:
        directory("pangenome/assemblies")
    run:
        for f in input:
            fout=os.path.basename(f)
            with gzip.open(f, 'rb') as f_in:
                with open(output[0] + "/" + fout[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        

rule aggregate_prokka:
    input_aggregate_forProkaa


# run prokkaC on downloaded genomes
# run roary on downloaded genomes

rule roary:
    input:
        "annotation/intermediary/prokka.faa.{db}.hmmscan",
        "annotation/annotation.filt.gff",
        "annotation/prokka.faa"
    output:
        "annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    params:
        indir="pangenome/annotation/",
        lambda wildcards: config["hmm_settings"][wildcards.db]["trim"]
    resources:
        runtime = "4:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_annotation.yaml"
    log: "logs/analysis_makegff.{db}.log"
    message: "makegff: Adding hmmer results of {wildcards.db} to gff."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        {SCRIPTSDIR}/hmmscan_addBest2gff.pl -i {input[0]} -a {input[1]} \
         -n {params[0]} {params[1]} -o {output} -g $(grep ">" {input[2]} | wc -l) > {log} 2>&1
        """

rule mergegff:
    input:
        "annotation/annotation.filt.gff",
        expand("annotation/intermediary/annotation.CDS.RNA.{db}.gff",db=config["hmm_DBs"].split())
    output:
        "annotation/annotation_CDS_RNA_hmms.gff"
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    threads: 1
    log: "logs/analysis_mergegff.log"
    message: "mergegff: Merging gffs with hmmer results."
    shell:
        """
        {SCRIPTSDIR}/mergegffs.pl {output} {input} > {log} 2>&1
        """

rule tar_annotation:
    input:
        "annotation/annotation_CDS_RNA_hmms.gff"
    output:
        "annotation/intermediary.tar.gz",
        "status/annotation.done"
    threads: 1
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    params:
        intermediary = "annotation/intermediary/"
    log: "logs/annotation_tar_annotation.log"
    message: "tar_annotation: Compressing intermediary steps of annotation."
    shell:
       """
       tar cvzf {output[0]} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
       touch {output[1]}
       """

