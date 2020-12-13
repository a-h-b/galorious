# set configuration file 
configfile: srcdir("config/config.default.yaml")

# set most important directories
SCRIPTSDIR = srcdir("workflow/scripts/")
ENVDIR = srcdir("workflow/envs/")
ROOTDIR = srcdir("/")

# include configuration file
include:
    "workflow/rules/get_config.smk"

# set working directory to output
workdir:
    OUTPUTDIR

# dump configuration into output directory
f = open('full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)


# useful functions
def getThreads(max):
    if workflow.cores:
        realThreads = max if max <= workflow.cores else workflow.cores
    elif workflow.nodes:
        realThreads = max if max <= workflow.nodes else workflow.nodes
    else:
        realThreads = max
    return realThreads

def filtering_input(wildcards):
    n = len(wildcards.filterstep.split("."))
    if n == 1:
        return ['Trimming/Illumina_fastq/r1.trimmed.fq',
        'Trimming/Illumina_fastq/r2.trimmed.fq',
        'Trimming/Illumina_fastq/se.trimmed.fq']
    elif n > 1:
        return [s + ".".join([s + "_filtered" for s in FILTER][:(n-1)]) + ".fq" for s in ['Trimming/Illumina_fastq/r1.trimmed.','Trimming/Illumina_fastq/r2.trimmed.','Trimming/Illumina_fastq/se.trimmed.' ]]
    else:
        raise ValueError("invalid filter length %s" % n)

def filtering_filter(wildcards):
    return ancient(expand(
            "{dir}/filtering/{filter}.{ext}", filter=wildcards.filterstep.split(".")[-1].split("_")[0],
            ext=['fa', 'fa.amb', 'fa.ann', 'fa.bwt', 'fa.pac', 'fa.sa'], dir=DBPATH))

# include rules based on chosen steps
final = []
if 'basecalling' in STEPS:
    include:
        "workflow/rules/long_read_basecalling.smk"
    final.append('status/basecalling.done')

if 'demultiplexing' in STEPS:
    include:
        "workflow/rules/long_read_demultiplexing.smk"
    final.append('status/demultiplexing.done')

if 'filtering_Nanopore' in STEPS:
    include:
        "workflow/rules/long_read_filtering.smk"
    final.append('status/filtering_Nanopore.done')

if 'filtering_Illumina' in STEPS:
    include:
        "workflow/rules/short_read_filtering.smk"
    final.append('status/filtering_Illumina.done')

if 'assembly' in STEPS:
    include:
        "workflow/rules/assembly.smk"
    final.append('status/assembly.done')

if 'annotation' in STEPS:
    include:
        "workflow/rules/annotation.smk"
    final.append('status/annotation.done')

if 'taxonomy_check' in STEPS:
    include:
        "workflow/rules/taxonomy.smk"
    final.append('status/taxonomy_check.done')

if 'pangenomics' in STEPS:
    include:
        "workflow/rules/pangenomics.smk"
    final.append('status/pangenomics.done')

# clean up in the end
if EMAIL == "":
    onsuccess:
        shell("mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv slurm* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv snakejob.* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log")
else:
    onsuccess:
        shell('mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv slurm* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; (  mv snakejob.* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; echo "$(date) {config[sessionName]}" | mail -s "galorious finish" {EMAIL} ')
    onerror:
        shell('echo "$(date) {config[sessionName]}" | mail -s "galorious exit with error" {EMAIL} ')
    onstart:
        shell('echo "$(date) {config[sessionName]}" | mail -s "galorious start" {EMAIL} ')


# master command
localrules: ALL
rule ALL:
    input:
        final
    output:
        touch("status/all.done") 





