# -*- coding: utf-8 -*-
"""
@created: Tue Oct 06 2020
@modified: Apr 27 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions for producing submatrices of the expression data matrices. Lighter for testing or even useful for some
analyses.
    _save_rna_tcga_test
    _generate_rna_tcga_test
    _save_rna_other_test
    _generate_rna_other_test

    generate_rna_test
"""

import os
import pandas as pd
import numpy as np

from pyprism.data import load_cln, load_rna_gex
from ._util_load import load_from_data, save_df_to_data
from ._filepaths import _get_filepath_res

DataFrame = pd.core.frame.DataFrame
Array = np.ndarray
Series = pd.core.series.Series

def _raise_unsupported_option(opt_user, opt_name, opt_poss):
    raise ValueError("Unsupported value %s for option '%s'. Choose one of %s" % (opt_user, opt_name, opt_poss))

def _random_indices(indices, n_rows: int, seed: int):
    np.random.seed(seed)
    idx_rand = np.random.permutation(len(indices))
    return indices[idx_rand[:n_rows]]

# TCGA =================================================================================================================

def _save_rna_tcga_test(df_test: DataFrame, metric: str):
    filepath = "./rna/tcga/tcga.Kallisto.%s_test.txt" % metric
    save_df_to_data(df_test, filepath, sep="\t", index=False)

def _generate_rna_tcga_test(rand_indices: int, n_cols: int, seed: int) -> DataFrame:
    # load
    df_rna_tcga_tpm = load_rna_gex(study="tcga", metric="tpm")
    df_rna_tcga_counts = load_rna_gex(study="tcga", metric="counts")

    # random indexing
    np.random.seed(seed)
    col_rand = np.random.permutation(df_rna_tcga_tpm.shape[1])
    columns = [x for x in df_rna_tcga_tpm.columns if x.startswith("TCGA")]
    rand_columns = columns[:n_cols]

    df_rna_tcga_tpm = df_rna_tcga_tpm.set_index("ensembl_gene_id")
    df_rna_tcga_counts = df_rna_tcga_counts.set_index("ensembl_gene_id")
    df_rna_tcga_tpm_test = df_rna_tcga_tpm.loc[rand_indices, rand_columns]
    df_rna_tcga_counts_test = df_rna_tcga_counts.loc[rand_indices, rand_columns]

    df_rna_tcga_tpm_test = df_rna_tcga_tpm_test.reset_index(drop=False)
    df_rna_tcga_tpm_test = df_rna_tcga_tpm_test.rename(columns={"ensembl_gene_id": "geneID"})
    df_rna_tcga_counts_test = df_rna_tcga_counts_test.reset_index(drop=False)
    df_rna_tcga_counts_test = df_rna_tcga_counts_test.rename(columns={"ensembl_gene_id": "geneID"})

    # save
    _save_rna_tcga_test(df_rna_tcga_tpm_test, "TPM")
    _save_rna_tcga_test(df_rna_tcga_counts_test, "counts")

# PRISM & MET500 =======================================================================================================

def _save_rna_other_test(df_test: DataFrame, study: str) -> None:
    filepath = "./rna/%s/quantification_test.csv" % study
    save_df_to_data(df_test, filepath, sep=",", index=False)

def _generate_rna_other_test(n_rows: int, study: str) -> DataFrame:
    # load
    df_rna = load_rna_gex(study=study, metric=None)
    df_rna = df_rna.set_index("ensembl_gene_id")

    # random indexing
    df_rna_test = df_rna.loc[rand_indices]
    df_rna_test = df_rna_test.reset_index(drop=False)

    # save
    _save_rna_other_test(df_rna_test, study)


# Function =============================================================================================================

def generate_rna_test(n_rows: int, n_cols_tcga: int, seed: int=0) -> None:
    indices = load_from_data(_get_filepath_res("gencode", "gencode_v27"))
    rand_indices = _random_indices(indices, n_rows, seed)

    _generate_rna_tcga_test(rand_indices, n_cols_tcga, seed)
    _generate_rna_other_test(rand_indices, "prism")
    _generate_rna_other_test(rand_indices, "met500")
