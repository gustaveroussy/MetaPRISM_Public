# -*- coding: utf-8 -*-
"""
@created: Apr 15 2021
@modified: Jan 06 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions for making metadata tables about sequencing data (wes and rna-seq).
"""

import os
import numpy  as     np
import pandas as     pd
import re
from typing import Iterable
import sys

from pyprism.data import load_ids

from ._filepaths import _get_filepath_rna_gex, _get_filepath_summary_rna_gex
from ._filepaths import _get_filepath_rna_fus, _get_filepath_summary_rna_fus
from ._filepaths import _get_filepath_wes_mut, _get_filepath_summary_wes_mut
from ._util_load import check_filepath_exists_in_data, load_from_data, save_df_to_data

#### type aliases
DataFrame = pd.core.frame.DataFrame

# RNA ==================================================================================================================

def _get_summary_rna_gex(study: str, level: str, metric: str, test_mode: str=False, other_mode: str=None) -> DataFrame:
    filepath = _get_filepath_rna_gex(study=study, level=level, metric=metric, test_mode=test_mode, other_mode=other_mode)
    cols_all = load_from_data(filepath, nrows=1).columns

    cols_gex = [x for x in cols_all if x not in ["geneID", "ensembl_gene_id"]]
    df_sum = pd.DataFrame({"Col_Name": cols_gex, "Sample_Id": cols_gex})
    df_ids = load_ids(study=study)
    return df_sum


def load_summary_rna_gex(study: str, level: str, metric: str, test_mode: bool=False, other_mode: str=None,
                     use_cache: bool=True) -> DataFrame:
    """
    Create a summary table linking Col_Name, Sample_Id and Subject_Id. The `Col_Name` is the raw name of the column in
    the expression table associated to the correspongind Sample_Id.

    You may specify a study as described below.

        - "met500": Trim Galore > Kallisto > Tximport. See https://github.com/gustaveroussy/MetaPRISM
        - "prism": Trim Galore > Kallisto > Tximport. See https://github.com/gustaveroussy/MetaPRISM
        - "tcga": From https://stanfordmedicine.app.box.com/s/lu703xuaulfz02vgd2lunxnvt4mfvo3q.
        - "tcga_6_samples_prism": PRISM pipeline on 6 TCGA FASTQs.
        - "tcga_6_samples_gao": Gao pipeline on 6 TCGA FASTQs.

    Parameters
    ----------
    See description in :func:`~pyprism.data.load_rna_gex`

    Returns
    -------
    df: DataFrame
    """
    filepath_cache = _get_filepath_summary_rna_gex(study, level, metric, test_mode, other_mode)
    filepath_cache_exists = check_filepath_exists_in_data(filepath_cache)

    if not filepath_cache_exists or not use_cache:
        # message
        print("-the cache file does not exist or use_cache was set to False, building it ... ", end="", flush=True)
        df_sum = _get_summary_rna_gex(study, level, metric, test_mode, other_mode)
        print("ok!", flush=True)
        if filepath_cache is not None:
            save_df_to_data(df_sum, filepath_cache, index=False)
    else:
        print("-using cache file for summary rna gex", flush=True)
        df_sum = load_from_data(filepath_cache)

    return df_sum



def _get_summary_rna_fus(study: str, mode: str=None) -> DataFrame:
    fields_id = ["Tumor_Sample_Barcode", "sample_ID", "Sample", "Sample_Id"]
    filepath = _get_filepath_rna_fus(study=study, mode=mode)
    cols_all = load_from_data(filepath, nrows=1).columns
    fields_id = list(set(cols_all).intersection(set(fields_id)))

    df_sum = load_from_data(filepath, usecols=fields_id)
    df_sum.columns = ["Sample_Id"]

    # get fusion counts
    df_sum = df_sum.groupby("Sample_Id").size().to_frame("N_Fusions").reset_index()

    # add subject ids
    df_ids = load_ids(study=study)
    df_sum = df_sum.merge(df_ids[["Sample_Id", "Subject_Id"]], how="left", on="Sample_Id")

    return df_sum


def load_summary_rna_fus(study: str, mode: str=None, use_cache: bool=True) -> DataFrame:
    """
    Create a summary table of the data in the aggregated fusion tables.

    Parameters
    ----------
    See description in :func:`~pyprism.data.load_rna_gex`

    use_cache: bool
        Set to False to force regeneration of the summary file.

    Returns
    -------
    df: DataFrame
    """
    filepath_cache = _get_filepath_summary_rna_fus(study, mode)
    filepath_cache_exists = check_filepath_exists_in_data(filepath_cache)

    if not filepath_cache_exists or not use_cache:
        # message
        print("-the cache file does not exist or use_cache was set to False, building it ... ", end="", flush=True)
        df_sum = _get_summary_rna_fus(study, mode)
        print("ok!", flush=True)
        if filepath_cache is not None:
            save_df_to_data(df_sum, filepath_cache, index=False)
    else:
        print("-using cache file for summary rna fus", flush=True)
        df_sum = load_from_data(filepath_cache)

    return df_sum

# WES ==================================================================================================================

def _get_summary_wes_mut(study: str, mode: str) -> DataFrame:
    if mode=="somatic_filters":
        fields_id_raw = ["Tumor_Sample", "Normal_Sample"]
        ids_old2new = {"Tumor_Sample": "Tumor_Sample_Barcode", "Normal_Sample": "Matched_Norm_Sample_Barcode"}
        fields_id = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    elif mode=="germline_maf":
        fields_id_raw = ["Sample_Id"]
        ids_old2new = None
        fields_id = ["Sample_Id"]
    else:
        fields_id_raw = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
        ids_old2new = None
        fields_id = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]


    # load only required columns
    filepath = _get_filepath_wes_mut(study, mode)
    df_sum = read_table(filepath, sep="\t", usecols=fields_id_raw)
    if ids_old2new is not None:
        df_sum = df_sum.rename(columns=ids_old2new)

    # get unique combinations of tumor - normal sample ids
    df_sum = df_sum.groupby(fields_id).size().to_frame("N_Mutations").reset_index()

    # add subject ids
    df_ids = load_ids(study=study)

    if "somatic" in mode:
        # add possible cases of samples without any detected event
        filepath = "%s/wes/somatic_maf/sample_list.tsv" % study
        df_sam = read_table(filepath)

        if study == "tcga":
            df_sam = df_sam.loc[df_sam["Comment"].fillna("NA").apply(lambda x: "MC3 Excluded" not in x)].copy()

        df_sum = df_sum.rename(columns={"Tumor_Sample_Barcode": "Tumor_Sample_Id",
                                        "Matched_Norm_Sample_Barcode": "Normal_Sample_Id"})
        cols = ["Tumor_Sample_Id", "Normal_Sample_Id"]
        df_sam = df_sam[cols]

        df_sum = df_sam.merge(df_sum, how="left")
        df_sum = df_sum.fillna(0)
        df_sum["N_Mutations"] = df_sum["N_Mutations"].astype(int)

        for nature in ["Tumor", "Normal"]:
            col_sample_id = "%s_Sample_Id" % nature

            df_sum = df_sum.merge(df_ids[["Sample_Id", "Subject_Id"]].drop_duplicates(),
                                  how="left",
                                  left_on=col_sample_id,
                                  right_on="Sample_Id")

            df_sum = df_sum.rename(columns = {"Subject_Id": "%s_Subject_Id" % nature})
            del df_sum["Sample_Id"]

        # check that subject id is identical between tumor and normal samples
        assert df_sum["Tumor_Subject_Id"].equals(df_sum["Normal_Subject_Id"])
        del df_sum["Tumor_Subject_Id"]
        df_sum = df_sum.rename(columns={"Normal_Subject_Id": "Subject_Id"})
    else:
        df_sum = df_sum.merge(df_ids[["Sample_Id", "Subject_Id"]], how="left", on="Sample_Id")

    return df_sum


def load_summary_wes_mut(study: str, mode: str, use_cache: bool=True):
    """
    Create a table linking Tumor_Sample_Id, Normal_Sample_Id and Subject_Id. It additionally computes the number of rows
    (i.e mutations) found in the mutation file for each pair of (Tumor_Sample_Id, Normal_Sample_Id).

    Parameters
    ----------
    See description in :func:`~pyprism.data.load_wes_mut`

    Returns
    -------
    df: A dataframe

    """
    filepath_cache = _get_filepath_summary_wes_mut(study, mode)
    filepath_cache_exists = check_filepath_exists_in_data(filepath_cache)

    if not filepath_cache_exists or not use_cache:
        # message
        print("-the cache file does not exist or use_cache was set to False, building it ... ", end="", flush=True)
        df_sum = _get_summary_wes_mut(study, mode, filepath_cache)
        print("ok!", flush=True)
        save_df_to_data(df_sum, filepath_cache, index=False)
    else:
        print("-using cache file for summary wes mut", flush=True)
        df_sum = load_from_data(filepath_cache)

    return df_sum
