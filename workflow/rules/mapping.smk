SAMTOOLS_MEM = str(round(float(config['bigMem'][:-1]) * 0.75 - 0.5)) + "G"

rule mapping_on_assembly:
    input:
        'trimming/Illumina_fastq/r1.preprocessed.fq',
        'trimming/Illumina_fastq/r2.preprocessed.fq',
        'trimming/Illumina_fastq/se.preprocessed.fq',
        "assembly/unicycler/assembly.fasta.amb",
        "assembly/unicycler/assembly.fasta.bwt",
        "assembly/unicycler/assembly.fasta.pac",
        "assembly/unicycler/assembly.fasta.sa",
        "assembly/unicycler/assembly.fasta.ann",
        "assembly/unicycler/assembly.fasta"
    output:
        'assembly/unicycler/reads.on.assembly.sorted.bam'
    resources:
        runtime = "12:00:00",
        mem = config['bigMem']
    params:
        prefix="assembly/unicycler/reads.on.assembly"
    threads: 2
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_mapping_on_assembly.log"
    message: "mapping_on_assembly: Mapping reads on merged assembly."
    shell:
        """
        SAMHEADER="@RG\\tID:{config[genomeName]}\\tSM:Illumina"
        PREFIX={params.prefix}
        USETHREADS="$(({threads}-1))"
        # merge paired and se
        samtools merge --threads $USETHREADS -f $PREFIX.merged.bam \
         <(bwa mem -v 1 -t $USETHREADS -M -R \"$SAMHEADER\" {input[8]} {input[0]} {input[1]} 2>> {log}| \
         samtools view --threads $USETHREADS -bS -) \
         <(bwa mem -v 1 -t $USETHREADS -M -R \"$SAMHEADER\" {input[8]} {input[2]} 2> {log}| \
         samtools view --threads $USETHREADS -bS -) 2>> {log}
        # sort
        samtools sort --threads $USETHREADS -m {SAMTOOLS_MEM} $PREFIX.merged.bam > $PREFIX.sorted.bam 2>> {log}
        rm $PREFIX.merged.bam
        """


rule bwa_index:
    input:
        "assembly/unicycler/assembly.fasta"
    output:
        "assembly/unicycler/assembly.fasta.amb",
        "assembly/unicycler/assembly.fasta.bwt",
        "assembly/unicycler/assembly.fasta.pac",
        "assembly/unicycler/assembly.fasta.sa",
        "assembly/unicycler/assembly.fasta.ann"
    resources:
        runtime = "24:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_bwa_index.log"
    message: "bwa_index: Indexing assembly for bwa."
    shell:
        """
        bwa index {wildcards.fasta} > {log} 2>&1
        """

rule index_bam:
    input:
        'assembly/unicycler/reads.on.assembly.sorted.bam'
    output:
        'assembly/unicycler/reads.on.assembly.sorted.bam.bai'
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_index.bam.log"
    message: "index_bam: Indexing bam."
    shell:
        """
        samtools index {input} > {log} 2>&1
        """


rule call_contig_depth:
    input:
        "assembly/unicycler/reads.on.assembly.sorted.bam",
        "assembly/unicycler/assembly.fasta",
        "assembly/unicycler/assembly.fasta.fai",
        "assembly/unicycler/assembly.fasta.bed3"
    output:
        "assembly/unicycler/assembly.contig_coverage.txt",
        "assembly/unicycler/assembly.contig_depth.txt",
        report("stats/assembly.contig_flagstat.txt",category="Assembly")
    resources:
        runtime = "2:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/analysis_call_contig_depth.log"
    message: "call_contig_depth: Getting data on assembly coverage with Illumina reads."
    shell:
        """
        coverageBed -b {input[0]} -a {input[3]} -sorted > {output[0]} 2>> {log}
        echo "Coverage calculation done" >> {log}
        echo "Running BEDTools for average depth in each position" >> {log}
        TMP_DEPTH=$(mktemp --tmpdir={TMPDIR} -t "depth_file_XXXXXX.txt")
        genomeCoverageBed -ibam {input[0]} | grep -v "genome" > $TMP_DEPTH
        echo "Depth calculation done" >> {log}

        ## This method of depth calculation was adapted and modified from the CONCOCT code
	perl {SCRIPTSDIR}/calcAvgCoverage.pl $TMP_DEPTH {input[1]} >{output[1]}	

        echo "Remove the temporary file" >> {log}
        rm $TMP_DEPTH
        echo "flagstat" >> {log}
        samtools flagstat {input[0]} 2>> {log} | cut -f1 -d ' ' > {output[2]}
        """

rule contig_fasta_indexing:
    input:
        "assembly/unicycler/assembly.fasta"
    output:
        "assembly/unicycler/assembly.fasta.fai"
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_contig_fasta_indexing.log"
    message: "contig_fasta_indexing: Indexing assembly."
    shell:
        """
        samtools faidx {input[0]} > {log} 2>&1
        """

rule contig_fasta2bed_conversion:
    input:
        "assembly/unicycler/assembly.fasta",
        "assembly/unicycler/assembly.fasta.fai"
    output:
        "assembly/unicycler/assembly.fasta.bed3"
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    threads: 1
    message: "contig_fasta2bed_conversion: Writing bed file for assembly."
    shell:
        """
        cat {input[0]}.fai | awk '{{print $1 \"\t0\t\" $2}}' > {output}
        """

