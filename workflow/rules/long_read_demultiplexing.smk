# run guppy
rule guppy_demulti:
    input:
        "basecalling/Nanopore_fastq"
    output:
        directory("demultiplexing/Nanopore_fastq/barcode%s" % config['barcode']),
        "status/demultiplexing.done"
    threads: 1
    log: "logs/guppy_demultiplexing.log"
    resources:
        runtime="24:00:00",
        mem=config['normalMem'],
    params:
        output="demultiplexing/Nanopore_fastq"
    shell:
        """
        {config[guppy][bin]}/guppy_barcoder --input_path {input} --save_path {params.output} --device auto --barcode_kits "{config[guppy][barcoding_kit]}" && touch {output[1]}
        """

