rule unicycler:
    input:
        "trimming/Nanopore_fastq/trimmed.fastq",
        expand('trimming/Illumina_fastq/{read}.preprocessed.fq',read=["r1","r2","se"]) 
    output:
        directory("assembly/unicycler"),
        "assembly/unicycler/assembly.fasta",
        "status/assembly.done"
    threads: 4
    log: "logs/assembly_unicycler.log"
    resources:
        runtime="24:00:00",
        mem=config['bigMem'],
    conda:
        ENVDIR + "galorious_assembly.yaml"
    params:
        keep=config['unicycler']['keep'],
        mode=config['unicycler']['mode'],
        linear=config['unicycler']['linear_seqs'],
        min=config['unicycler']['min_fasta_length']
    shell:
        """
        unicycler -1 {input[1]} -2 {input[2]} -s {input[3]} -l {input[0]} -o {output[0]} --keep {params.keep} --mode {params.mode} --linear_seqs {params.linear} --min_fasta_length {params.min} -t {threads} --verbosity 3 &>> {log} && touch {output[2]} 
        """
