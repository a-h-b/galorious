SAMTOOLS_MEM = str(round(float(config['bigMem'][:-1]) * 0.75 - 0.5)) + "G"

FILTERING_CMD = """
        TMP_FILE=$(mktemp --tmpdir={TMPDIR} -t "alignment_XXXXXX.bam")
        bwa mem -v 1 -t {threads} {input[3]} {input[0]} {input[1]} 2>> {log} | samtools view --threads {threads} -bS - > $TMP_FILE 2> {log}
        samtools merge --threads {threads} -u - \
         <(samtools view --threads {threads} -u  -f 4 -F 264 $TMP_FILE 2>> {log}) \
         <(samtools view --threads {threads} -u -f 8 -F 260 $TMP_FILE 2>> {log}) \
         <(samtools view --threads {threads} -u -f 12 -F 256 $TMP_FILE 2>> {log}) 2>> {log} | \
         samtools view --threads {threads} -bF 0x800 - 2> {log} | \
         samtools sort --threads {threads} -m {SAMTOOLS_MEM} -n - 2> {log} | \
         bamToFastq -i stdin -fq {output[0]} -fq2 {output[1]} >> {log} 2>&1
        if [[ -s {input[2]} ]]
         then
         bwa mem -v 1 -t {threads} {input[3]} {input[2]} 2>> {log} | \
          samtools view --threads {threads} -bS - 2>> {log} | \
          samtools view --threads {threads} -uf 4 - 2>> {log} | \
          bamToFastq -i stdin -fq {output[2]} >> {log} 2>&1
        else
         echo "{input[2]} is empty, skipping single end reference sequence filtering, but creating it anyway..." >> {log}
         touch {output[2]}
        fi
        rm -rf $TMP_FILE
        """

rule illumina_filtering:
    input:
        filtering_input,
        filtering_filter
    output:
        expand('trimming/Illumina_fastq/{read}.trimmed.{{filterstep}}.fq',read=["r1","r2","se"])
    wildcard_constraints:
        filterstep = "|".join([".".join([s + "_filtered" for s in FILTER][:x]) for x in range(1,len(FILTER)+1)])
    resources:
        runtime = "24:00:00",
        mem = config['bigMem']
    threads: getThreads(config['bigCores'])
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/filtering.illumina.{filterstep}.log"
    message: "filtering: Filtering Illumina reads to get {wildcards.filterstep}."
    shell:
        FILTERING_CMD


rule symlink_illumina_preprocessed_files:
    input:
        expand("trimming/Illumina_fastq/{read}.trimmed.{filtered}.fq",filtered=".".join([s + "_filtered" for s in FILTER]),read=["r1","r2","se"])
    output:
        expand('trimming/Illumina_fastq/{read}.preprocessed.fq',read=["r1","r2","se"])
    resources:
        runtime = "1:00:00",
        mem = config['normMem']
    threads: 1
    message: "symlink_trimmed_filtered_files: Symlinking filtered and trimmed reads."
    run:
        for i,o in zip(input,output):
             shell("ln -fs $(echo {i} | cut -f 2 -d /) {o} && touch -h {o}")


if config['nextseq']:
    rule trimming:
        input:
            'trimming/Illumina_fastq/r1.fq',
            'trimming/Illumina_fastq/r2.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'trimming/Illumina_fastq/r1.trimmoed.fq',
            'trimming/Illumina_fastq/se1.trimmoed.fq',
            'trimming/Illumina_fastq/r2.trimmoed.fq',
            'trimming/Illumina_fastq/se2.trimmoed.fq'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = config['normMem']
        conda: ENVDIR + "galorious_trimming.yaml"
        log: "logs/trimming.illumina.log"
        message: "trimming illumina reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """

    rule nextseq_trimming:
        input:
            'trimming/Illumina_fastq/{read}.trimmoed.fq'
        output:
            'trimming/Illumina_fastq/{read}.trimmed.fq',
            'trimming/Illumina_fastq/{read}.trimmoed.fq.gz'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = config['normMem']
        conda: ENVDIR + "galorious_trimming.yaml"
        log: "logs/trimming.illumina.nextseq.{read}.log"
        message: "trimming: Trimmming poly-G's off {wildcards.read} reads."
        shell:
            """
            cutadapt -j {threads} --nextseq-trim={config[trimmomatic][trailing]} -o {output[0]} {input} >> {log} 2>&1
            pigz -p {threads} {input} >> {log} 2>&1
            """
else:
    rule trimming:
        input:
            'trimming/Illumina_fastq/r1.fq',
            'trimming/Illumina_fastq/r2.fq',
            DBPATH + "/adapters/adapters.done"
        output:
            'trimming/Illumina_fastq/r1.trimmed.fq',
            'trimming/Illumina_fastq/se1.trimmed.fq',
            'trimming/Illumina_fastq/r2.trimmed.fq',
            'trimming/Illumina_fastq/se2.trimmed.fq'
        threads: getThreads(10)
        resources:
            runtime="12:00:00",
            mem = config['normMem']
        conda: ENVDIR + "galorious_trimming.yaml"
        log: "logs/trimming.illumina.log"
        message: "trimming: Trimmming illumina reads."
        shell:
            """
            trimmomatic PE -threads {threads} {input[0]} {input[1]} {output} \
             ILLUMINACLIP:{DBPATH}/adapters/{config[trimmomatic][adapter]}.fa:{config[trimmomatic][seed_mismatch]}:{config[trimmomatic][palindrome_clip_threshold]}:{config[trimmomatic][simple_clip_threshold]} \
             LEADING:{config[trimmomatic][leading]} \
             TRAILING:{config[trimmomatic][trailing]} \
             SLIDINGWINDOW:{config[trimmomatic][window_size]}:{config[trimmomatic][window_quality]} \
             MINLEN:{config[trimmomatic][minlen]} \
             MAXINFO:{config[trimmomatic][target_length]}:{config[trimmomatic][strictness]} > {log} 2>&1
            """


rule cat_se_trimmed:
    input:
        'trimming/Illumina_fastq/se1.trimmed.fq',
        'trimming/Illumina_fastq/se2.trimmed.fq'
    output:
        'trimming/Illumina_fastq/se.trimmed.fq'
    threads: 1
    resources:
        runtime="12:00:00",
        mem = config['normMem']
    message: "cat_se_trimmed: Concatenating trimmed single end reads"
    shell:
        """
        cat {input[0]} {input[1]} > {output}
        """

localrules: ctrl_short_trim
rule ctrl_short_trim:
    input:
        expand('trimming/Illumina_fastq/{read}.preprocessed.fq',read=["r1","r2","se"])
    output:
        touch('status/filtering_Illumina.done')

