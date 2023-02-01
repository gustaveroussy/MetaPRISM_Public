# -*- coding: utf-8 -*-
"""
@author: Yoann Pradat

Tests for _utils.py module
"""

import pandas as pd
from pyprism.util import explode_df

def test_explode_df():
    df_test = pd.DataFrame({"A": [0,1,2], "B": ["a,b", "c", "d"]})
    df_expl = explode_df(df_test, cols="B", sep=",", fill_value="")
    df_good = pd.DataFrame({"A": [0,0,1,2], "B": ["a", "b", "c", "d"]})

    assert df_expl.equals(df_good)
