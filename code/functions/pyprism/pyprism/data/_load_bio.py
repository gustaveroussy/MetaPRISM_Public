# -*- coding: utf-8 -*-
"""
@created: Nov 27 2020
@modified: Dec 15 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions for loading biospecimen annotations files.
    1. One file with biospecimens annotations for TCGA
    2. One file with biospecimens annotations for PRISM

Functions
    _load_bio_prism

    load_bio
"""

import os
import pandas as pd
from typing import Iterable

from ._util_load import load_from_data, subset_data
from ._filepaths import _get_filepath_bio_curated

DataFrame = pd.core.frame.DataFrame

# User function ========================================================================================================

def load_bio(study: str, identifiers: Iterable[str]=None, identifiers_name: str=None, mode: str="in_design", **kwargs):
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
        Choose "in_design" to load only data for samples in the design or 'all' to load data for all samples of the
        cohort.

    Returns
    -------
    df_bio: DataFrame
    """
    df_bio = load_from_data(_get_filepath_bio_curated(study, mode), low_memory=False, **kwargs)
    df_bio = subset_data(df_bio, values=identifiers, col_name=identifiers_name)

    return df_bio
