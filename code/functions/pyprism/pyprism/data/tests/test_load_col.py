# -*- coding: utf-8 -*-
"""
@created: May 19 2021
@modified: May 19 2021
@author: Yoann Pradat

Tests for _load_colors.py module.
"""

import pandas as pd
from pyprism.data import load_colors

# User functions =======================================================================================================

def test_load_colors():
    global2colors = load_colors(sheet="Biopsy_Site")
    assert type(global2colors)==dict

    df_global_colors = load_colors(sheet="Biopsy_Site", as_dataframe=True)
    assert type(df_global_colors) == pd.core.frame.DataFrame

    BS2colors = load_colors(sheet="Biopsy_Site")
    assert type(BS2colors)==dict




