# -*- coding: utf-8 -*-
"""
@created: Apr 22 2021
@created: Oct 01 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Code for the user function load_cln.
"""

import os
import numpy as np
import pandas as pd
import re
import sys

from typing import Iterable
from pyprism.util import setwd_to_scripts

cwd = setwd_to_scripts()
sys.path.append("pipeline_cln/workflow/functions")
from cols_2_groups import _get_cols_2_groups_cln_brca_prism, _get_cols_2_groups_cln_brca_tcga
from cols_2_groups import _get_cols_2_groups_cln_met500, _get_cols_2_groups_cln_prism, _get_cols_2_groups_cln_tcga
os.chdir(cwd)

from ._util_load import load_from_data, subset_data, make_column_groups
from ._filepaths import _get_filepath_cln_curated

DataFrame = pd.core.frame.DataFrame

def _get_cols_2_groups(study: str, mode: str):
    if study == "met500":
        return _get_cols_2_groups_cln_met500(mode)
    elif study == "prism":
        if mode in ["all", "in_design"]:
            return _get_cols_2_groups_cln_prism(mode)
        elif mode=="brca":
            return _get_cols_2_groups_cln_brca_prism()
        else:
            raise ValueError("Unsupported combination of 'study' and 'mode'")
    elif study == "tcga":
        if mode in ["all", "in_design"]:
            return _get_cols_2_groups_cln_tcga(mode)
        elif mode=="brca":
            return _get_cols_2_groups_cln_brca_tcga()
        else:
            raise ValueError("Unsupported combination of 'study' and 'mode'")


# User function ========================================================================================================

def load_cln(study: str, identifiers: Iterable[str]=None, identifiers_name: str=None, multi_cols=False,
             mode:str="in_design", **kwargs) -> DataFrame:
    """
    Load clinical tables for the specified study.

    Parameters
    ----------
    study: str
        Choose one of 'tcga' or 'prism'.
    identifiers: Iterable[str], optional.
        If not None, return annotations only for the identifiers values.
    identifiers_name: Iterable[str], optional.
        Name of the column of identifiers.
    multi_cols: bool, optional
        Choose False to load only the last level of column names.
    mode: str, optional
        Either 'in_design' or 'all' for 'pancancer' annotations or 'brca' for per-cancer annotations.

    Returns
    -------
    df_cln: DataFrame
    """

    df_cln = load_from_data(_get_filepath_cln_curated(study=study, mode=mode), **kwargs)
    df_cln = subset_data(df_cln, values=identifiers, col_name=identifiers_name)

    if multi_cols:
        cols2groups = _get_cols_2_groups(study, mode)
        df_cln = make_column_groups(df_cln, cols2groups)

    return df_cln
