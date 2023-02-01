# -*- coding: utf-8 -*-
"""
@created: Oct 06 2020
@modified: Aug 30 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions for loading in one-line different expression data matrices.
    load_rna_fus(study: str, mode: str=None, identifiers=None, identifiers_name=None):
    load_rna_gex(study: str, level: str=None, metric: str="counts", test_mode: bool=False, other_mode: str=None,
             identifiers=None, identifiers_name=None):
"""

from   functools import reduce
import traceback
import os
import numpy  as     np
import pandas as     pd
from   typing import Dict, Iterable
import re

from pyprism.data import load_summary_rna_gex
from ._filepaths import _get_filepath_rna_gex, _get_filepath_rna_fus
from ._util_load import load_from_data, subset_data

#### type aliases
DataFrame = pd.core.frame.DataFrame

def load_rna_fus(study: str, mode: str=None, identifiers: Iterable[str]=None, identifiers_name: str=None) -> DataFrame:
    """
    You may specify a study and a mode as described below.

        - "met500": arriba, ericscript, pizzly, starfusion, oncokb_civic
        - "prism": arriba, ericscript, pizzly, starfusion, oncokb_civic
        - "tcga": aggregated (Gao supp table), oncokb_civic

    Parameters
    ----------
    study:  str
        Name of the cohort used in the naming of the files. Choose one of: "met500", "prism" or "tcga".
    model: str
        You may specify 'arriba', 'ericscript', 'pizzly', 'starfusion' or 'oncokb_civic' for 'met500' or 'prism'. For
        'tcga', only 'aggregated' (from Gao et al paper) and 'oncokb_civic' are valid entries.
    identifiers: Iterable[str], optional.
        If not None, return RNA data only for the identifiers values.
    identifiers_name: str, optional
        Field name of the identifiers.

    Returns
    -------
    df: DataFrame
        The data loaded in a DataFrame
    """
    # load maf file
    print("-reading aggregated RNA table ... ", end="")
    filepath = _get_filepath_rna_fus(study, mode)
    df_rna = load_from_data(filepath)
    print("ok!")

    df_rna = subset_data(df_rna, values=identifiers, col_name=identifiers_name)
    return df_rna



def load_rna_gex(study: str, level: str="genes", metric: str="counts", test_mode: bool=False, other_mode: str=None,
                 identifiers: Iterable[str]=None, identifiers_name: str=None) -> DataFrame:
    """
    You may specify a study as described below.

        - "met500": Trim Galore > Kallisto > Tximport. See https://github.com/gustaveroussy/MetaPRISM
        - "prism": Trim Galore > Kallisto > Tximport. See https://github.com/gustaveroussy/MetaPRISM
        - "tcga": From https://stanfordmedicine.app.box.com/s/lu703xuaulfz02vgd2lunxnvt4mfvo3q.
        - "tcga_6_samples": PRISM/Gao pipeline on 6 TCGA FASTQs.

    Parameters
    ----------
    study:  str
        Name of the cohort used in the naming of the files. Choose one of: "prism", "tcga".
    level:  str, optional
        Level at which counts are aggregated. You may choose "genes" or "transcripts".
    metric: str, optional
        Choose "TPM", "counts" or "length". Default is 'counts'.
    test_mode: bool, optional
        If set to True, load submatrices generated with generate_gex_test from util_rna.py module.
    other_mode: str, optional
        Use 'prism' or 'zheng' when study='tcga_6_samples'.
    identifiers: Iterable[str], optional.
        If not None, return RNA data only for the identifiers values.
    identifiers_name: str, optional
        Field name of the identifiers.

    Returns
    -------
    df: DataFrame
        The data loaded in a DataFrame
    """
    # load metadata
    df_summary = load_summary_rna_gex(study, level=level, metric=metric, test_mode=test_mode, other_mode=other_mode)

    # select on identifiers
    if identifiers is not None and identifiers_name is not None:
        df_summary = df_summary.loc[df_summary[identifiers_name].isin(identifiers)].copy()

    usecols = df_summary["Col_Name"].tolist()

    # load gex file
    print("-reading aggregated RNA table ... ", end="")
    filepath = _get_filepath_rna_gex(study, level, metric, test_mode, other_mode)
    if study=="tcga":
        if level=="genes":
            usecols = ["geneID"] + usecols
        df_rna = load_from_data(filepath, usecols=usecols)
        df_rna = df_rna.rename(columns={"geneID": "ensembl_gene_id"})
    else:
        if level=="genes":
            usecols = ["ensembl_gene_id"] + usecols
        else:
            usecols = ["target_id"] + usecols
        df_rna = load_from_data(filepath, usecols=usecols)
    print("ok!")

    return df_rna
