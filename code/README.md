# Code organisation

## functions

This folder hosts the code for R and python packages that were developed for this project as well as links to external
tools. In short,

- `rprism` R package with useful functions for reading, writing, and manipulating data tables of the project as well as
  performing tasks repeated in many analyses such as computation of pvalues and multiple testing correction.

- `rprismtools` R package with useful functions for drawing figures. 

- `pyprism` Python package with useful functions for reading, writing, and manipulating data tables of the project. It
  provides python implementations for many of the functions available in `rprism`.

- `pyprismtools` Python package with useful functions for drawing figures.

- `tools` See README files of the linked submodules.

**Contributors**
- [Yoann Pradat](https://github.com/ypradat)

# pipelines

This folder hosts the bioinformatics pipelines that were developed or edit from existing pipelines in order to run all
bioinformatics analyses for the project. In short,

- `kallisto_tximport` a Snakemake port of [RNASeq_pipeline](https://github.com/gevaertlab/RNASeq_pipeline) by
  [gevaertlab](https://github.com/gevaertlab). Command line were duplicated as-is, additional logging files were added
  to follow the pipeline's execution cautiously. See `pipelines/kallisto_tximport/README.md` for more details.

- `rnafusion` clone of the [nf-core rnafusion pipeline](https://github.com/nf-core/rnafusion) v1.2.0.
  See `pipelines/rnafusion/README.md` for more details.

- `wes` a Snakemake pipeline developed specifically for the project. It performs a comprehensive analysis of WES data in
  order to identify, filter, and annotate germline mutations, somatic mutations, copy-number alterations, purity and
  ploidy, and microsatellite instability in each pair of tumor/normal. It may also run in tumor-only mode for almost
  all analyses except microsatellite instability. See `pipelines/wes/README.md` for more details.


**Contributors**
- [Thibault Dayris](https://github.com/tdayris)
- [Marc Deloger](https://github.com/mdeloger)
- [Ismael Padioleau](https://github.com/ipadiole)
- [Yoann Pradat](https://github.com/ypradat)

# scripts

This folder gather all the scripts that were used to perform the analysis, extract the numbers, and draw the figures
that underlie the paper. Each analysis is organised in a subfolder using a structure similar to that of the
`pipelines/wes` pipeline. The code is organized and made easily runnable through Snakemake. See `scripts/README.md` for
more details.

**Contributors**
- [Daniel Gautheret](https://github.com/dgautheret)
- Konstantin Gunbin
- [Antoine Lainé](https://github.com/aLaine1)
- [Yoann Pradat](https://github.com/ypradat)
- Andrey Yurchenko
