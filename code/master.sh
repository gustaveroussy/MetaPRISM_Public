#!/usr/bin/env bash

# activate snakemake
export PATH="${HOME}/miniconda3/bin:${HOME}/miniconda3/condabin:$PATH"
source activate snakemake

# echo current path
path=`pwd`
printf "the current path is: %s\n" "$path"

# move to an analysis folder - "data_overview"
cd scripts/data_overview

file_conda_log=../common/logs/setup_conda.done
if [[ -f "${file_conda_log}" ]]; then
  rm ${file_conda_log}
fi

# snakemake command for reproducing F1a_left.svg
snakemake -s workflow/Snakefile \
  --use-conda \
  --conda-prefix ${HOME}/.config/snakemake/conda \
  --jobs 1 ../../../results/figures_paper/F1a_left.svg -f
