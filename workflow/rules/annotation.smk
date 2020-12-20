rule prokkaC:
    input:
        "assembly/unicycler/assembly.fasta" 
    output:
        "annotation/annotation.filt.gff",
        "annotation/prokka.gff",
        "annotation/prokka.faa",
        "annotation/prokka.fna",
        "annotation/prokka.ffn",
        "annotation/prokka.fsa"
    threads: getThreads(8)
    log: "logs/analysis_annotate.log"
    resources:
        runtime="4:00:00",
        mem=config['normalMem']
    params:
        out="annotation",
        pref="prokka"
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

        # Prokka gives a gff file with a long header and with all the contigs at the bottom.  The command below removes the
        # And keeps only the gff table.

        LN=`grep -Hn "^>" {output[1]} | head -n1 | cut -f2 -d ":" || if [[ $? -eq 141 ]]; then true; else exit $?; fi`
        LN1=1
        LN=$(($LN-$LN1))
        head -n $LN {output[1]} | grep -v "^#" | sort | uniq | grep -v "^==" > {output[0]}
        """

rule hmmer:
    input:
        "annotation/prokka.faa"
    output:
        "annotation/intermediary/prokka.faa.{db}.hmmscan"
    params: 
        lambda wildcards: config["hmm_settings"][wildcards.db]["cutoff"],
        dbs = DBPATH + "/hmm/{db}"
    resources:
        runtime = "2:00:00",
        mem = config['normalMem'] 
    conda: ENVDIR + "galorious_annotation.yaml"
    threads: getThreads(12)
    log: "logs/analysis_hmmer.{db}.log"
    message: "hmmer: Running HMMER for {wildcards.db}."
    shell:
        """
         hmmsearch --cpu {threads} {params[0]} --noali --notextw \
          --tblout {output} {params.dbs}/*.hmm {input} >/dev/null 2> {log}
        """

rule makegff:
    input:
        "annotation/intermediary/prokka.faa.{db}.hmmscan",
        "annotation/annotation.filt.gff",
        "annotation/prokka.faa"
    output:
        "annotation/intermediary/annotation.CDS.RNA.{db}.gff"
    params:
        "{db}",
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

