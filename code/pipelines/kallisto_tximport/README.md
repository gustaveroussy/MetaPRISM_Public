# MetaPRISM pipelines

This repository holds pipelines for MetaPRISM project. Please find the following repositories:

-   `envs`: contains environment and tools versions
-   `pipelines`: contains snakefiles and other pipeline scripts
-   `reports`: contains descriptions of outputs provided by each pipeline
-   `scripts`: contains additional scripts that were developed by for the purpose of this project
-   `tests`: contains small test data-sets in order to guarantee the well-being of the pipelines


## Pipeline Kallisto-tximport

This pipeline is a Snakemake port of [RNASeq_pipeline](https://github.com/gevaertlab/RNASeq_pipeline) by [gevaertlab](https://github.com/gevaertlab). Command line were duplicated as-is, additional logging files were added to follow the pipeline's execution cautiously.

To run Kallisto-tximport, please refer to [snakemake](https://snakemake.readthedocs.io/en/stable/index.html)'s documentation. In order to run this pipeline on the IGR's cluster, please refer to our dedicated [snakemake profile](https://github.com/tdayris-perso/slurm).

If you have any question, please contact MetaPRISM team or IGR's bioinformatics core platform.

Warning: NONCODEv5 annotation is not available. Use the GencodeV27 annotation.

Please note that:

-   Annotations are downloaded and prepared by this pipeline on-the-fly
-   Tools (with versions provided in the [original pipeline](https://github.com/gevaertlab/RNASeq_pipeline)) are installed by this pipeline on-the-fly
-   Singularity container is provided with this pipeline, pulled and built on-the-fly
-   Cold/hot storage as defined in Flamingo cluster at the IGR is taken into account by this pipeline, no need to copy fastq files by yourself. Let the pipeline do it automatically. These copy are automatically deleted once no more needed.


An example of command line for local run would be:

      snakemake -s /path/to/Kallisto_tximport.smk --use-conda --use-singularity --reason --printshellcmds

An exemple of command line on IGR's Flamingo would be:

      snakemake -s /path/to/Kallisto_tximport.smk --profile slurm
