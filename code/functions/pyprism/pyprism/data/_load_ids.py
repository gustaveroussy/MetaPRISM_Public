# -*- coding: utf-8 -*-
"""
@created: 04/15/21
@modified: 04/15/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Load the curated tables of identifiers.
"""

import os
import pandas as pd

from ._util_load import load_from_data
from ._filepaths import _get_filepath_ids_curated

DataFrame = pd.core.frame.DataFrame

# User function ========================================================================================================

def load_ids(study: str):
    """
    Load clinical tables for the specified study.

    Parameters
    ----------
    study: str
        Choose one of 'met500, 'prism' or 'tcga'.

    Returns
    -------
    df_ids DataFrame
    """
    return load_from_data(_get_filepath_ids_curated(study), dtype="str")
