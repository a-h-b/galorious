localrules: ctrl_anno
def ctrl_anno_input(steplist):
    input = ["annotation/checkM/checkM.result.tsv",
             "annotation/intermediary.tar.gz"]
    if "assembly" in steplist:
        input.append(["stats/assembly.contig_flagstat.txt",
                      "assembly/reads.on.assembly.sorted.bam.bai"])
    return input

rule ctrl_anno:
    input:
        ctrl_anno_input(STEPS)
    output:
        touch("status/annotation.done")


rule prokkaC:
    input:
        "assembly/assembly.fasta"
    output:
        "annotation/annotation.filt.gff",
        "annotation/prokka.gff",
        "annotation/prokka.faa",
        "annotation/prokka.fna",
        "annotation/prokka.ffn",
        "annotation/prokka.fsa"
    threads: getThreads(4)
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
    threads: getThreads(4)
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
        runtime = "2:00:00",
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
        "annotation/annotation_antiSmash.gff",
        expand("annotation/intermediary/annotation.CDS.RNA.{db}.gff",db=config["hmm_DBs"].split())
    output:
        "annotation/annotation_CDS_RNA_hmms.gff"
    resources:
        runtime = "2:00:00",
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
        "annotation/intermediary.tar.gz"
    threads: 1
    resources:
        runtime = "2:00:00",
        mem = config['normalMem']
    params:
        intermediary = "annotation/intermediary/"
    log: "logs/annotation_tar_annotation.log"
    message: "tar_annotation: Compressing intermediary steps of annotation."
    shell:
        """
        tar cvzf {output[0]} {params.intermediary} >> {log} 2>&1 && rm -r {params.intermediary} >> {log} 2>&1
        """


if "taxonomy_check" in STEPS:
    rule checkM:
        input:
            "taxonomy/GTDB/gtdbtk.bac120.summary.tsv",
            "annotation/prokka.faa"
        output:
            "annotation/checkM/checkM.result.tsv"
        threads: 1
        resources:
            runtime = "2:00:00",
            mem = config['normalMem']
        params:
            annodir = "annotation",
            checkmdir = "annotation/checkM"
        message: "checkM: from prokkaC annotation for GTDB species"
        log: "logs/checkM.log"
        conda: ENVDIR + "galorious_check.yaml"
        shell:
            """
            TAX=`cut -f 2 {input[0]} | tail -n 1`
            mkdir -p {params.checkmdir}
            checkm taxon_list >> {params.checkmdir}/taxa.txt
            IFS=';' read -ra TAXA <<< "$TAX"
            IFS=',' read -ra RANKS <<< "domain,phylum,class,order,family,genus,species"
            for i in `seq $((${{#TAXA[@]}}-1)) -1 0`
            do
               RANK=${{RANKS[$i]}}
               TAXON=`echo ${{TAXA[i]}} | sed "s#^.__##" | sed "s#_.##g"`
               if grep "$RANK" -w {params.checkmdir}/taxa.txt | grep "$TAXON" -w -q; then 
                   checkm taxon_set "$RANK" "$TAXON" {params.checkmdir}/${{TAXON// /}}.ms &>> {log}
                   checkm analyze -x faa -g {params.checkmdir}/${{TAXON// /}}.ms \
                    {params.annodir} {params.checkmdir} &>> {log}
                   checkm qa -o 2 -f {output} --tab_table {params.checkmdir}/${{TAXON// /}}.ms \
                    {params.checkmdir} &>> {log}
                   break
               fi
            done
            """
else:
    rule checkM_lineage:
        input:
            "annotation/prokka.faa"
        output:
            "annotation/checkM/checkM.result.tsv"
        threads: 2
        resources:
            runtime = "2:00:00",
            mem = config['bigMem']
        params:
            annodir = "annotation",
            checkmdir = "annotation/checkM"
        message: "checkM: from prokkaC annotation"
        log: "logs/checkM.log"
        conda: ENVDIR + "galorious_check.yaml"
        shell:
            """
            checkm lineage_wf --genes -t {threads} --tab_table -f {output} -x faa {params.annodir} {params.checkmdir} &>> {log}
            """ 

checkpoint antismash:
    input:
        "assembly/assembly.fasta",
        "annotation/annotation.filt.gff"
    output:
        directory("annotation/antiSmash")
    params:
        db = DBPATH + "/antismash_db"
    resources:
        runtime = "4:00:00",
        mem = config['normalMem']
    threads: getThreads(4)
    log: "logs/annotation_antismash.log"
    message: "antismash: Running antiSmash on the assembly."
    conda:
        ENVDIR + "galorious_antiSmash.yaml"
    shell:
        """
        antismash -c {threads} --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --databases {params.db} --genefinding-gff3 {input[1]} --genefinding-tool none --output-dir {output} {input[0]} &>> {log}
        """

rule gbk2gff:
    input:
        "annotation/antiSmash/{somefile}.gbk"
    output:
        "annotation/antiSmash/{somefile}.gbk.gff"
    threads: 1
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    params:
        outdir="annotation/antiSmash"
    log: "logs/gbk2gff.{somefile}.log"
    message: "gbk2gff: converting {wildcards.somefile}"
    conda: 
        ENVDIR + "galorious_roary.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        perl bp_genbank2gff3.pl --outdir {params.outdir} {input} &>> {log}
        IN={input}
        keep=$(basename ${{IN%%.*}})
        sed -i "/^[^$keep]/d" {output}
        """

def gather_regions(wildcards):
        checkpoint_output=checkpoints.antismash.get().output[0]
        all = expand("annotation/antiSmash/{i}.gbk.gff",
                                 i=glob_wildcards(os.path.join(checkpoint_output,"{i}.gbk")).i)
        all.remove("annotation/antiSmash/assembly.gff")
        return all

rule regions2gff:
    input:
        gff="annotation/annotation.filt.gff",
        regions=gather_regions
    output:
        "annotation/annotation_antiSmash.gff"
    params:
        gffs = lambda wildcards,input: input['regions'].join(",")
    threads: 1
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    log: "logs/regions2gff.log"
    message: "regions2gff: adding antiSmash regions to gff"
    conda:
        ENVDIR + "galorious_annotation.yaml"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        export LC_ALL=en_US.utf-8
        perl {SCRIPTSDIR}/addAntiSmash2gff.pl -a {input[0]} -s {params.gffs} -o {output}
        """

