# -*- coding: utf-8 -*-
"""
@created: Apr 15 2021
@modified: Apr 15 2021
@author: Yoann Pradat

Tests for _load_ids.py module.
"""

import pandas as pd
from .._load_ids import load_ids

# User functions =======================================================================================================

def test_load_ids():
    df_ids = load_ids(study="met500")
    assert type(df_ids) == pd.core.frame.DataFrame

    df_ids = load_ids(study="prism")
    assert type(df_ids) == pd.core.frame.DataFrame

    df_ids = load_ids(study="tcga")
    assert type(df_ids) == pd.core.frame.DataFrame
