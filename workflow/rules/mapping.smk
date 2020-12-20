SAMTOOLS_MEM = str(round(float(config['normalMem'][:-1]) * 0.75 - 0.5)) + "G"

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
        mem = config['normalMem']
    params:
        prefix="assembly/unicycler/reads.on.assembly"
    threads: getThreads(4)
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_mapping_on_assembly.log"
    message: "mapping_on_assembly: Mapping reads on merged assembly."
    shell:
        """
        SAMHEADER="@RG\\tID:{config[genomeName]}\\tSM:Illumina"
        PREFIX={params.prefix}
        # merge paired and se
        samtools merge --threads {threads} -f $PREFIX.merged.bam \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[0]} {input[1]} 2>> {log}| \
         samtools view --threads {threads} -bS -) \
         <(bwa mem -v 1 -t {threads} -M -R \"$SAMHEADER\" {input[8]} {input[2]} 2> {log}| \
         samtools view --threads {threads} -bS -) 2>> {log}
        # sort
        samtools sort --threads {threads} -m {SAMTOOLS_MEM} $PREFIX.merged.bam > $PREFIX.sorted.bam 2>> {log}
        rm $PREFIX.merged.bam
        """


rule bwa_index:
    input:
        "{fasta}"
    output:
        "{fasta}.amb",
        "{fasta}.bwt",
        "{fasta}.pac",
        "{fasta}.sa",
        "{fasta}.ann"
    resources:
        runtime = "24:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "/IMP_mapping.yaml"
    log: "logs/assembly_bwa_index.{fasta}.log"
    message: "bwa_index: Indexing {wildcards.fasta} for bwa."
    shell:
        """
        bwa index {wildcards.fasta} > {log} 2>&1
        """

rule index_bam:
    input:
        '{sorted_bam}'
    output:
        '{sorted_bam}.bai'
    resources:
        runtime = "8:00:00",
        mem = config['normalMem']
    threads: 1
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/assembly_index.{sorted_bam}.log"
    message: "index_bam: Indexing {wildcards.sorted_bam}."
    shell:
        """
        samtools index {input} > {log} 2>&1
        """



