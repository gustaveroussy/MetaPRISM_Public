# -*- coding: utf-8 -*-
"""
@created: 05/19/21
@modified: 16/06/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Function for loading lists of colors.
"""

import os
import numpy as np
import pandas as pd

from pyprism.util import setwd_to_data
from ._filepaths import _get_filepath_colors

DataFrame = pd.core.frame.DataFrame


# User function ========================================================================================================

def load_colors(sheet: str, as_dataframe=False):
    """Load one or multiple sheets from an excel workbook containing lists of colors for categories of different fields.

    Parameters
    ----------
    sheet: str
        Name of the excel sheet to be loaded.
    as_dataframe: bool, optional
        Set to True if you would like to have the list of colors returned in a dataframe. If False, colors are returned
        as a dictionary

    Returns
    -------
    A dict or a data.frame.

    """
    try:
        filepath = _get_filepath_colors()
        cwd = setwd_to_data()
        df_colors = pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl", keep_default_na=False)
        os.chdir(cwd)
    except Exception as e:
        raise e

    if as_dataframe:
        return df_colors
    else:
        colors = dict(zip(df_colors["Name"], df_colors["Color"]))
        return colors
