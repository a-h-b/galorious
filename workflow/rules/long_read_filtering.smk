rule nanofilt:
    input:
        "demultiplexing/Nanopore_fastq/barcode%s" % config['barcode']
    output:
        "trimming/Nanopore_fastq/trimmed.fastq",
        "status/filtering_Nanopore.done"
    threads: 1
    log: "logs/filtering_Nanopore.log"
    resources:
        runtime="24:00:00",
        mem=config['normalMem'],
    conda:
        ENVDIR + "galorious_trimming.yaml"
    params:
        qcutoff=config['nanofilt']['quality'],
        length=config['nanofilt']['length']
    shell:
        """
        cat {input}/*.fastq | Nanofilt -q {params.qcutoff} -l {params.length} > {output[0]} && touch {output[1]}
        """
