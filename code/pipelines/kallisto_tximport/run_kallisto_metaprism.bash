#!/usr/bin/bash
source /mnt/beegfs/software/miniconda/4.8.3/etc/profile.d/conda.sh
conda activate /mnt/beegfs/software/snakemake/5.20.1
module load singularity
snakemake -s /mnt/beegfs/scratch/m_deloger/MetaPRISM/RNAseq/pipeline/MetaPRISM/pipelines/Kallisto_tximport/Kallisto_tximport.smk --use-conda --use-singularity --reason --printshellcmds --profile slurm --configfile /mnt/beegfs/scratch/m_deloger/MetaPRISM/RNAseq/pipeline/MetaPRISM/pipelines/Kallisto_tximport/config.gencode.yaml --keep-going
