# -*- coding: utf-8 -*-
"""
@modified: Dec 22 2020
@created: Dec 22 2020
@author: Yocln Pradat

Tests for _load_cln.py module.
"""

import pandas as pd
from pyprism.data import LoaderCln

def test_LoaderCln():
    def _date_birth_keep(x):
        return pd.to_datetime(x["Date_Birth"]) > pd.to_datetime("1960-01-01")

    loader = LoaderCln(
        studies=["prism"],
        cln2keep_each = {
            "prism": {"Project_TCGA": "LUAD"},
            "tcga": {"Project_TCGA": "LUAD"},
        },
        cln2apply_all = {
            ("Date_Birth",): ("age below 60",  _date_birth_keep)
        },
        verbose=True)

    df_cln = loader.load()
