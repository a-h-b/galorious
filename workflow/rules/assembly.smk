def input_unicycler_longreads(step_list, input):
    if "filtering_Nanopore" in step_list:
        return "trimming/Nanopore_fastq/trimmed.fastq"
    elif "basecalling" in step_list or "demultiplexing" in step_list:
        raise Exception("You've skipped filtering_Nanopore but did basecalling and/or demultiplexing. This is not recommended. If you don't want to loose data in the filtering step, set quality and length to low values in the config, but do include filtering_Nanopore in steps.")
    else:
        return input

def input_unicycler_shortreads(step_list, input):
    if "filtering_Illumina" in step_list:
        return expand('trimming/Illumina_fastq/{read}.preprocessed.fq',read=["r1","r2","se"])
    else:
        return input

rule unicycler:
    input:
        input_unicycler_longreads(STEPS,config['inputs']['Nanopore']),
        input_unicycler_shortreads(STEPS,config['inputs']['Illumina'])
    output:
        "assembly/unicycler/assembly.fasta",
        "assembly/assembly.fasta"
    threads: 4
    log: "logs/assembly_unicycler.log"
    resources:
        runtime="24:00:00",
        mem=config['bigMem'],
    conda:
        ENVDIR + "galorious_assembly.yaml"
    params:
        outdir="assembly/unicycler",
        keep=config['unicycler']['keep'],
        mode=config['unicycler']['mode'],
        linear=config['unicycler']['linear_seqs'],
        min=config['unicycler']['min_fasta_length']
    shell:
        """
        unicycler -1 {input[1]} -2 {input[2]} -s {input[3]} -l {input[0]} \
         -o {params.outdir} --keep {params.keep} --mode {params.mode} --linear_seqs \
         {params.linear} --min_fasta_length {params.min} -t {threads} --verbosity 3 &>> {log}
        sed 's/ .*//g' {output[0]} > {output[1]} 
        """

localrules: assembly_ctrl

rule assembly_ctrl:
    input:
        "stats/assembly.contig_flagstat.txt",
        "assembly/reads.on.assembly.sorted.bam.bai"
    output:
        touch("status/assembly.done")


