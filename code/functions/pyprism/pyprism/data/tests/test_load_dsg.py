# -*- coding: utf-8 -*-
"""
@created: Feb 12 2021
@modified: Feb 12 2021
@author: Yoann Pradat

Tests for _load_dsg.py module.
"""

import pandas as pd
from .._load_dsg import load_dsg

# User functions =======================================================================================================

def test_load_dsg():
    df_design = load_dsg(study="prism")
    assert type(df_design) == pd.core.frame.DataFrame
