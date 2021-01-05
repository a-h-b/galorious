localrules: cp_assembly_as_input
rule cp_assembly_as_input:
    input:
        config['inputs']['Contigs']
    output:
        "assembly/assembly.fasta"
    run:
        shutil.copy(input, output)
