# Rprism

Rprism is an R package that supports the analysis performed in the META-PRISM project. The following sections provide 
installation instructions and examples of how to use the package.

## Installation

### Build

If the package is not available as a compressed tarball (`.tar.gz` file), you need to build this file before you can 
install the package. You may build it by running the command `make build`.

### Install

Use the latest tarball available at the root to install the package. In case you are using `renv` to handle R libraries
within your project, you may copy this tarball in the `renv/local` folder and then simply run the command
`install.packages("rprism")`. This command will automatically install all dependencies listed in the `DESCRIPTION` file
and the current package. Please make sure that the `install.packages` command can find all dependencies by specifying
the repositories (CRAN, BioCsoft, BiocCann, BioCexp, etc.) where it should look for packages. Once again, if you use
`renv` you can specify default repositories in the `renv.lock` file.

## Examples

R package that contains the following
- `load_xxx.R`: exports functions to load data.
    Examples

    ```
    library(rprism)

    # Samples data
    help(load_bio)

    df_bio <- load_bio(study="prism")
    df_bio <- load_bio(study="tcga")

    # Clinical data
    help(load_cln)

    df_cln <- load_cln(study="prism")
    df_cln <- load_cln(study="tcga", identifiers=c("TCGA-02-0047", "TCGA-ZS-A9CG"), identifiers_name="Subject_Id")
    df_cln <- load_cln(study="tcga", mode="brca")

    # Fusions data
    help(load_fus)
    df_fus <- load_fus(study="met500", mode="unfiltered")
    df_fus <- load_fus(study="prism", mode="filtered")
    df_fus <- load_fus(study="tcga", mode="filtered")

    # RNA data (this may take a few mintues to run)
    help(load_rna_gex)
    df_rna <- load_rna_gex(study="met500", metric="counts")
    df_rna <- load_rna_gex(study="prism", mode="TPM")
    df_rna <- load_rna_gex(study="tcga", mode="length")

    # WES data (this may take a few mintues to run)
    help(load_wes_mut)

    df_maf <- load_wes_mut(study="met500", mode="somatic_filtered")
    df_maf <- load_wes_mut(study="met500", mode="somatic_pathogenic")
    df_maf <- load_wes_mut(study="prism", mode="germline_unfiltered")
    df_maf <- load_wes_mut(study="tcga", mode="somatic_filtered")

    # WES and RNA metadata
    help(load_metadata_rna)
    df_meta <- load_metadata_rna(study="met500")
    df_meta <- load_metadata_rna(study="prism")
    df_meta <- load_metadata_rna(study="tcga")

    help(load_metadata_wes)
    df_meta <- load_metadata_wes(study="met500", mode="somatic_filtered")
    df_meta <- load_metadata_wes(study="prism", mode="germline_pathogenic")
    df_meta <- load_metadata_wes(study="tcga", mode="somatic_filtered")
    ```

### Testing

In order to run the tests of rprism, run

```
make test
```

Look into the documentation of 
[devtools::test()](https://www.rdocumentation.org/packages/devtools/versions/2.3.2/topics/test) for more details.
