# run guppy
rule guppy_demulti:
    input:
        "basecalling/Nanopore_fastq"
    output:
        directory("demultiplexing/Nanopore_fastq")
    threads: 1
    log: "logs/guppy_demultiplexing.log"
    resources:
        runtime="24:00:00",
        mem=config['normalMem'],
    shell:
        """
        {config[guppy][bin]}/guppy_barcoder --input_path {input} --save_path {output} --device auto --barcode_kits {config[guppy][barcoding_kit]}
        """

