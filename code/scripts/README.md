# Organisation

The code for each separate part of the analysis is located in the folders listed below along with a brief description of
the analysis performed in each folder.

- `combined_alterations`
    - Contains scripts for preparing tables for META-PRISM, MET500, and TCGA that combine annotated point mutations and
      small indels, annotated focal SCNAs, annotated fusions, annotated MSI status, annotated TMB status, and annotated
      AR-V7 expression score.
    - Hosts the code for drawing Fig. 3A, Fig. 3B, Fig. 5A, Fig. 5B, Supplementary Fig. S9, and Supplementary Fig. S10.
    - Contains scripts for preparing count tables of the number of DNA (WES-derived only) or DNA & RNA
    (WES- and RNAseq-derived) events per gene, per pathway, or in total.
- `common`
    - Contains scripts and files for setting up conda environments that were used for running the analyses.
    - Hosts the code to perform the selection of analyzable samples in each of META-PRISM, MET500, and TCGA for each
    type of analysis.
    - Hosts the scripts for the preparation and processing of aggregated tables with one line per patient and
      clinical and genomics variables in columns for use in modeling analyses.
    - Hosts the scripts for the plots detailing the variable selection process, the model qualities and the model
      coefficients.
- `data_overview`
    - Hosts the code for drawing Fig. 1A, Fig. 1B, Fig. 1C, Fig. 1D, Fig. 2B (middle), Fig. 5B, Supplementary Fig.
      S1, Supplementary Fig. S2, Supplementary Fig. S3A, Supplementary Fig. S3B, Supplementary Fig. S3C, and
      Supplementary Fig. S3D.
- `fusions_analysis`
    - Hosts the code that was used to perform the filtering of fusions and identify the best combination of algorithms
    in META-PRISM, MET500, and TCGA.
    - Hosts the code for drawing Fig. 4B and Supplementary Methods. Fig 2, Supplementary Methods Fig. 5, Supplementary
      Methods Fig. 10, and Supplementary Methods Fig. 11.
- `germline_mutations`
    - Hosts the code for drawing Fig. 3C, Fig. 3D and Fig. 3E.
- `immuno_analysis`
    - Contains scripts for computing ssGSEA scores of the 29 signatures in META-PRISM, MET500, and TCGA; training
      predictive models on TCGA scores; applying the models on META-PRISM, MET500 and TCGA in order to classify each
      RNAseq sample in one of the four tumor types from Bagaev *et al*. 2019.
    - Hosts the code for drawing Fig. 4A.
- `mutational_signatures`
    - Contains scripts for performing the mutational signature analysis on WES samples from META-PRISM, MET500, and
      TCGA.
    - Hosts the code for drawing Fig. 2B (top), Fig. 2B (bot), and Supplementary Fig. S4.
- `numbers_paper`
    - Contains scripts for computing all the numbers written in the text.
- `somatic_cnas`
    - Contains scripts for preparing count tables of the number of SCNA events per gene or per pathway.
    - Hosts the code for drawing Fig. 2C, Fig. 2D, Supplementary Fig. S6, and Supplementary Fig. S7.
- `somatic_mutations`
    - Contains scripts for preparing inputs to Mutpanning and for running it.
    - Hosts the code for drawing Fig. 2A and Supplementary Fig. S8.
    - Contains scripts for preparing count tables of the number of mutations per gene, per pathway, or in total.
- `survival_analysis`
    - Contains scripts for running survival models.
    - Hosts the code for drawing all Kaplan-Meier curves including Fig. 2E and Fig. 7A.
    - Hosts the code for Fig. 7C, Fig. 7D, Fig. 7E, Fig. 7F.
    - Hosts the code for drawing Supplementary Fig. S11.
- `treatment_resistances`
    - Contains scripts for processing manual annotations of associations between treatments and oncogenic alterations.
    - Hosts the code for drawing Fig. 6.

# Per analysis

## Organisation

The workflow of each analysis is divided in rules that are assembled in the snakemake subfiles of `workflow/rules` which
are themselves included in the global snakemake file `workflow/Snakefile`. The rules call scripts located in the
`workflow/scripts` folder.

The parameters of the analysis may be modified in the `config/config.yaml` file. All libraries needed for running the
analysis are installed during the first steps of the pipeline through [conda virtual
environments](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) whose specifications
are given in the folder `common/envs`.

## How to run the analysis?

In order to run the analysis, ensure the `snakemake` and `conda` commands are available (if you are on a HPC using
slurm, you may use `module load` to load [modules](https://curc.readthedocs.io/en/latest/compute/modules.html) providing
you with these commands). You can launch the full pipeline via

```
snakemake -s workflow/Snakefile --jobs [n_jobs] --profile [your_profile]
```

where `[n_jobs]` is the maximum number of CPUs (or jobs if you are on a cluster) to be used/run in parallel and
`[your_profile]` your snakemake profile (read the [snakemake
documentation](<https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles>)).

In case you are only interested in running part of the analysis so as to reproduce one of the results file, you may do
so with

```
snakemake -s workflow/Snakefile --jobs [n_jobs] --profile [your_profile] [path/to/file/you/want/to/reproduce]
```
