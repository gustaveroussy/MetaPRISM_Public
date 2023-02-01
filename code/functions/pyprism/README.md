# pyprism

## Installation

From an activated Python environment (with conda or any other environment manager), install the python package with 
(cd into the folder `functions/pyprism`). 

```
make install
```

### Description

Python package that contains the following subpackages

- `data`: functions for loading and performing formatting operations on data matrices.
    Examples

    ```
    from pyprism.data import load_cln
    from pyprism.data import load_bio
    from pyprism.data import load_design
    from pyprism.data import load_maf, load_summary_maf
    from pyprism.data import load_rna, load_summary_rna

    from pyprism.data import delineate_data

    help(load_cln)
    help(load_bio)
    help(load_design)
    help(load_maf)
    help(load_rna)

    # pan-cancer annotations
    df_cln = load_cln(study="met500")
    df_cln = load_cln(study="prism")
    df_cln = load_cln(study="tcga", identifiers=["TCGA-02-0047", "TCGA-ZS-A9CG"], identifiers_name="Subject_Id")

    # per-cancer annotations
    df_cln = load_cln(study="tcga", mode="brdca")

    # biospecimen data
    df_bio = load_bio(study="met500")
    df_bio = load_bio(study="prism")
    df_bio = load_bio(study="tcga")

    # WES (maf) data (this may take a few minutes to run)
    df_maf = load_maf(study="met500")
    df_maf = load_maf(study="prism")
    df_maf = load_maf(study="tcga")

    # summary
    df_meta = load_summary_rna(study="met500")
    df_meta = load_summary_rna(study="prism")
    df_meta = load_summary_rna(study="tcga")

    df_meta = load_summary_wes(study="met500", mode="somatic_unfiltered")
    df_meta = load_summary_wes(study="prism", mode="germline_pathogenic")
    df_meta = load_summary_wes(study="tcga", mode="somatic_filtered")

    # test_mode=True loads a subset of the data; loading the full table is memory-consuming for tcga.
    df_rna = load_rna(study="met500", metric="counts", test_mode=True)
    df_rna = load_rna(study="prism", metric="counts", test_mode=True)
    df_rna = load_rna(study="prism", metric="TPM", test_mode=True)
    df_rna = load_rna(study="tcga", metric="TPM", test_mode=True)
    ```

- `stats`: functions for running statistical methods.
    An extensive example is implemented in `tests/test_neighborhood.py`.

    ```
    from pyprism.stats import NeighbordhoodAnalysis

    neigh = NeighborhoodAnalysis(corr_type="binary", grid_lims=(-0.5, 1.5), grid_size=250, n_perm=5, seed=123, n_jobs=1)
    neigh.fit(X,y)
    ```

- `util`: functions for performing technical operations (e.g setting the working directory or saving results).
    Examples

    ```
    import os
    from pyprism.util import setwd_to_data, explode_df

    cwd = setwd_to_data()
    print("cwd before %s" % cwd)
    print("cwd after %s" % os.getcwd())
    ```

## Testing

Use the following command from the `functions/pyprism` folder to run the tests of pyprism

```
make test
```

In order to run this command successfully, you will need the data folder as it is on the nextcloud. Look into the
documentation of [pytest](https://docs.pytest.org/en/stable/) for more advanced testing options.
