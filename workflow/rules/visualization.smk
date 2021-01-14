
def circos_input(steplist):
    inputs = ["assembly/assembly.length.txt",
              "annotation/annotation_CDS_RNA_hmms.gff",
              "assembly/assembly.gc_skew.txt",
              "assembly/assembly.contig_depth_perBase.txt"]
    if "analyse_pangenome" in steplist:
        inputs.append("pangenome/roary/gene_presence_absence.anno.RDS")
    if "additional_mapping" in steplist:
        inputs.append("status/additional_mapping.done")
    return inputs

rule prep_circos:
    input: 
        circos_input(STEPS)
    output:
        "visualization/data/contigs.ideogram"
    threads: 1
    conda: ENVDIR + "galorious_visualization.yaml"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    log: "logs/prep_circos.log"
    script:
        SCRIPTSDIR + "prep_circos.R"

rule circos:
    input:
        "visualization/data/contigs.ideogram"
    output:
        "visualization/galorious.svg"
    threads: 1
    conda: ENVDIR + "galorious_visualization.yaml"
    resources:
        runtime="1:00:00",
        mem=config['normalMem']
    log: "logs/circos.log"
    shell:
        """
        export PERL5LIB=$CONDA_PREFIX/lib/site_perl/5.26.2
        circos -conf {SCRIPTSDIR}circos.conf &>> {log}
        """

