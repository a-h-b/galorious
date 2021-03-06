SAMTOOLS_MEM = str(round(float(config['bigMem'][:-1]) * 0.75 - 0.5)) + "G"

rule mapping_on_assembly:
    input:
        'trimming/Illumina_fastq/r1.preprocessed.fq',
        'trimming/Illumina_fastq/r2.preprocessed.fq',
        'trimming/Illumina_fastq/se.preprocessed.fq',
        "assembly/assembly.fasta.amb",
        "assembly/assembly.fasta.bwt",
        "assembly/assembly.fasta.pac",
        "assembly/assembly.fasta.sa",
        "assembly/assembly.fasta.ann",
        "assembly/assembly.fasta"
    output:
        'assembly/reads.on.assembly.sorted.bam'
    resources:
        runtime = "12:00:00",
        mem = config['bigMem']
    params:
        prefix="assembly/reads.on.assembly"
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
        ancient("assembly/assembly.fasta")
    output:
        "assembly/assembly.fasta.amb",
        "assembly/assembly.fasta.bwt",
        "assembly/assembly.fasta.pac",
        "assembly/assembly.fasta.sa",
        "assembly/assembly.fasta.ann"
    resources:
        runtime = "24:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_bwa_index.log"
    message: "bwa_index: Indexing assembly for bwa."
    shell:
        """
        bwa index {input} > {log} 2>&1
        """

rule index_bam:
    input:
        '{bam}.sorted.bam'
    output:
        '{bam}.sorted.bam.bai'
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_index.{bam}.log"
    message: "index_bam: Indexing {wildcards.bam}."
    shell:
        """
        samtools index {input} > {log} 2>&1
        """

rule contig_length:
    input:
        ancient("assembly/assembly.fasta")
    output:
        "assembly/assembly.length.txt",
        "assembly/assembly.gc_content.txt",
        "assembly/assembly.gc_skew.txt"
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_annotation.yaml"
    log: "logs/contig_length.log"
    message: "contig_length: Getting data on assembly length and GC content."
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        perl {SCRIPTSDIR}fastaNamesSizes.pl {input} > {output[0]} 2>> {log}

        echo "Obtaining GC content"
        TMP_GC=$(mktemp --tmpdir={TMPDIR} -t "gc_out_XXXXXX.txt")
        perl {SCRIPTSDIR}get_GC_content.pl {input} $TMP_GC >> {log} 2>&1

        # The program above provides a file gc_out.txt. This command cleans the output
        echo "Clean up output" >> {log}
        cut -f1,2 $TMP_GC | sed -e 's/>//g'> {output[1]}
        echo "Remove intermediate files" >> {log}
        rm $TMP_GC

        echo "GCskew"
        python {SCRIPTSDIR}gcskew.py -i {input} -o {output[2]} -k 2000 >>{log} 2>&1 
        """


rule call_contig_depth:
    input:
        "assembly/reads.on.assembly.sorted.bam",
        "assembly/assembly.fasta",
        "assembly/assembly.fasta.fai",
        "assembly/assembly.fasta.bed3",
        "assembly/assembly.length.txt"
    output:
        "assembly/assembly.contig_coverage.txt",
        "assembly/assembly.contig_depth.txt",
        "assembly/assembly.contig_depth_perBase.txt"
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
        genomeCoverageBed -ibam {input[0]} 2>> {log} | grep -v "genome" > $TMP_DEPTH 2>> {log}
        genomeCoverageBed -ibam {input[0]} -d -g {input[4]} >{output[2]} 2>> {log}  
        echo "Depth calculation done" >> {log}

        ## This method of depth calculation was adapted and modified from the CONCOCT code
	perl {SCRIPTSDIR}calcAvgCoverage.pl $TMP_DEPTH {input[1]} >{output[1]}	

        echo "Remove the temporary file" >> {log}
        rm $TMP_DEPTH
        """

rule contig_fasta_indexing:
    input:
        ancient("assembly/assembly.fasta")
    output:
        "assembly/assembly.fasta.fai"
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
        "assembly/assembly.fasta.fai"
    output:
        "assembly/assembly.fasta.bed3"
    resources:
        runtime = "1:00:00",
        mem = config['normalMem']
    threads: 1
    message: "contig_fasta2bed_conversion: Writing bed file for assembly."
    shell:
        """
        cat {input[0]} | awk '{{print $1 \"\t0\t\" $2}}' > {output}
        """

