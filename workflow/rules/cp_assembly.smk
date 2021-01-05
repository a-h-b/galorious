localrules: cp_assembly_as_input
rule cp_assembly_as_input:
    input:
        config['inputs']['Contigs'],
        config['inputs']['Illumina'].split(" ")
    output:
        "assembly/assembly.fasta",
        "trimming/Illumina_fastq/r1.preprocessed.fq",
        "trimming/Illumina_fastq/r2.preprocessed.fq",
        "trimming/Illumina_fastq/se.preprocessed.fq"
    run:
        for i, o in zip(input,output):
            shutil.copy(i, o)
