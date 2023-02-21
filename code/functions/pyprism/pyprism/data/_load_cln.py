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
from ._util_load import load_from_data, subset_data
from ._filepaths import _get_filepath_cln_curated

DataFrame = pd.core.frame.DataFrame

# User function ========================================================================================================

def load_cln(study: str, identifiers: Iterable[str]=None, identifiers_name: str=None, mode:str="in_design", **kwargs) -> DataFrame:
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
    mode: str, optional
        Either 'in_design' or 'all' for 'pancancer' annotations or 'brca' for per-cancer annotations.

    Returns
    -------
    df_cln: DataFrame
    """

    df_cln = load_from_data(_get_filepath_cln_curated(study=study, mode=mode), **kwargs)
    df_cln = subset_data(df_cln, values=identifiers, col_name=identifiers_name)

    return df_cln
