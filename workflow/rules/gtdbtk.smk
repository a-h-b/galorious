rule GTDBtk:
    input:
        "assembly/unicycler/assembly.fasta"
    output:
        "taxonomy/GTDB/gtdbtk.bac120.summary.tsv",
        "taxonomy/species",
        "status/taxonomy.done"
    params:
        dir="taxonomy/",
        outdir="taxonomy/GTDB"
    resources:
        runtime = "3:00:00",
        mem = config['bigMem']
    conda: ENVDIR + "galorious_gtdbtk.yaml"
    log: "logs/taxonomy_GTDBtk.log"
    message: "GTDBtk: Classifiying using GTDBtk."
    threads: getThreads(config['bigCores'])
    shell:
        """
        mkdir -p {params.dir} 
        cp {input} {params.dir}
        export GTDBTK_DATA_PATH="{DBPATH}/GTDB_tk"
        export PYTHONPATH=$CONDA_PREFIX/lib/python3.7/site-packages
        gtdbtk classify_wf --genome_dir {params.dir} -x fasta --out_dir {params.outdir} --cpus {threads} > {log} 2>&1
        cut -f 2 {output[0]} | tail -n 1  > {output[1]} 
        touch {output[2]}
        """

