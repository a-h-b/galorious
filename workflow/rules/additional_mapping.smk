import pandas as pd

if os.path.isabs(os.path.expandvars(config['inputs']['Additional_Illumina_samples'])):
    sam_path = os.path.expandvars(config['inputs']['Additional_Illumina_samples'])
else:
    sam_path = os.getcwd() + "/" + os.path.expandvars(config['inputs']['Additional_Illumina_samples'])
try:
    samples = pd.read_table(sam_path)
except:
    print("Table for additional Illumina samples was not found. Please enter the absolute path and file name in the config file or remove additional_mapping from steps list.")
    raise
if 'r1_file' not in samples.columns:
    raise Exception("You haven't provided file names for read 1 - column should be named r1_file.")
if 'r2_file' not in samples.columns:
    raise Exception("You haven't provided file names for read 2 - column should be named r2_file.")
if 'se_file' not in samples.columns:
    raise Exception("You haven't provided file names for single end reads - column should be named se_file.")
if 'path' not in samples.columns:
    raise Exception("You haven't provided a path to look for the additional samples - column should be named path.")
samples = samples.set_index("sample",drop=False)
samples.index = samples.index.set_levels([i.astype(str) for i in samples.index.levels])


def getAdditionalMapping_input(wildcards):
    return samples.loc[wildcards.sample, "path"] + "/" + samples.loc[wildcards.sample,["r1_file", "r2_file","se_file"].dropna()


rule add_mapping_on_assembly:
    input:
        getAdditionalMapping_input,
        "assembly/assembly.fasta.amb",
        "assembly/assembly.fasta.bwt",
        "assembly/assembly.fasta.pac",
        "assembly/assembly.fasta.sa",
        "assembly/assembly.fasta.ann",
        "assembly/assembly.fasta"
    output:
        'additional_mapping/{sample}.reads.on.assembly.sorted.bam'
    resources:
        runtime = "12:00:00",
        mem = config['bigMem']
    params:
        prefix="additional_mapping/{wildcards.sample}.reads.on.assembly"
    threads: 2
    conda: ENVDIR + "galorious_mapping.yaml"
    log: "logs/additional_mapping_on_assembly.{sample}.log"
    message: "mapping_on_assembly: Mapping {wildcards.sample} reads on merged assembly."
    shell:
        """
        SAMHEADER="@RG\\tID:{config[genomeName]}\\tSM:{wildcards.sample}"
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


rule call_contig_add_depth:
    input:
        "additional_mapping/{sample}.reads.on.assembly.sorted.bam",
        "assembly/assembly.fasta",
        "assembly/assembly.fasta.fai",
        "assembly/assembly.fasta.bed3"
    output:
        "additional_mapping/{sample}.contig_coverage.txt",
        "additional_mapping/{sample}.contig_depth.txt",
        report("additional_mapping/{sample}.contig_flagstat.txt",category="Assembly")
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


localrules: ctrl_addMap
rule ctrl_addMap
    input:
        expand("additional_mapping/{sam}.contig_depth.txt",sam=samples.sample)
    output:
        touch("status/additional_mapping.done")
