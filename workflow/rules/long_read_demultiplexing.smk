def input_guppy_demulti(step_list, input):
    if "basecalling" in step_list:
        return "basecalling/Nanopore_fastq"
    else:
        return input


# run guppy
rule guppy_demulti:
    input:
        input_guppy_demulti(STEPS,config['inputs']['Nanopore'])
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
    message: "demultiplexing"
    shell:
        """
        {config[guppy][bin]}/guppy_barcoder --input_path {input} --save_path {params.output} --device auto --barcode_kits "{config[guppy][barcoding_kit]}" && touch {output[1]}
        """

