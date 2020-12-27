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
* GTDB database
You can get this [as explained here](https://ecogenomics.github.io/GTDBTk/installing/index.html) like so:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
tar xvzf gtdbtk_r95_data.tar.gz
```
Give the path to the extracted directory in the config file.

You also need the taxonomy data from GTDB, which is in the taxonomy directory of the extracted directory.
If you don't want to set it in the config file and it will automatically download itself into your output directory.
In addition, you need the full metadata set from GTDB:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/bac120_metadata_r95.tar.gz
tar xvzf bac120_metadata_r95.tar.gz
```
If you don't want to download it, you can keep the field empty in the config file and it will automatically download itself.

~~* ncbi's datasets command line tool: it will install itself, if necessary.~~ download is done by entrez direct tools, installed via conda.
