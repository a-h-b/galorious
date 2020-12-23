rule GTDBtk:
    input:
        "assembly/unicycler/assembly.fasta"
    output:
        directory("taxonomy/GTDB"),
        "status/taxonomy.done"
    params:
        dir="taxonomy/"
    resources:
        runtime = "12:00:00",
        mem = config['bigMem']
    conda: ENVDIR + "galorious_gtdbtk.yaml"
    log: "logs/taxonomy_GTDBtk.log"
    message: "GTDBtk: Classifiying using GTDBtk."
    threads: 4
    shell:
        """
        mkdir -p {params.dir} 
        cp {input} {params.dir}
        export GTDBTK_DATA_PATH="{DBPATH}/GTDB_tk"
        export PYTHONPATH=$CONDA_PREFIX/lib/python3.7/site-packages
        gtdbtk classify_wf --genome_dir {params.dir} -x fasta --out_dir {output[0]} --cpus {threads} > {log} 2>&1
        touch {output[1]}
        """

