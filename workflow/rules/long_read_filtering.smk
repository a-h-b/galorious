def input_nanofilt(step_list, input, barcode):
    if "demultiplexing" in step_list:
        return "demultiplexing/Nanopore_fastq/barcode%s" % barcode
    elif "basecalling" in step_list:
        return "basecalling/Nanopore_fastq"
    else:
        return input



rule nanofilt:
    input:
        input_nanofile(STEPS,config['inputs']['Nanopore'],config['barcode'])
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
        cat {input}/*.fastq | NanoFilt -q {params.qcutoff} -l {params.length} > {output[0]} && touch {output[1]}
        """
