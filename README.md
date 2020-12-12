# galorious

galorious' purpose is to assemble and analyse bacterial isolate genomes from Nanopore and Illumina reads.

### Dependencies:
* conda
* snakemake
If you want to use snakemake via conda, install the environment, as [recommended by Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html):
```
conda install -c conda-forge mamba
mamba create --prefix $PWD/conda/snakemake_env
conda activate $PWD/conda/snakemake_env
mamba install -c conda-forge -c bioconda snakemake
conda deactivate
```
* ONT's guppy basecaller and barcoder
Available from Oxford Nanopore Technologies. If additional dependencies are necessary (e.g. CUDA), you can prepend them in the command. galorious expects GPU's to be used when submitting the basecalling and demultiplexing steps to a cluster.
