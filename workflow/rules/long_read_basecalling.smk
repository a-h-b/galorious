# run guppy
if config['guppy']['model']:
    GUPPY_SHELL: """
        {config[guppy][bin]}/guppy_basecaller --input_path {input} --save_path {output[0]} --device auto --config {config[guppy][model]} && touch {output[1]} 
    """
elif config['guppy']['flowcell']:
    GUPPY_SHELL: """
        {config[guppy][bin]}/guppy_basecaller --input_path {input} --save_path {output[0]} --device auto --flowcell {config[guppy][flowcell] --kit config[guppy][sequencing_kit]} && touch {output[1]}
    """

rule guppy_basecalling:
    input: 
        config['inputs']['Nanopore']
    output: 
        directory("basecalling/Nanopore_fastq"),
        "status/basecalling.done"
    threads: 1
    log: "logs/guppy_basecalling.log"
    resources:
        runtime="24:00:00",
        mem=config['normalMem'],
    shell:
        GUPPY_SHELL


# report stats on how many reads were in the fast5 files and how much was written?
# run fastqc?


