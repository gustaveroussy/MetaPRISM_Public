# -*- coding: utf-8 -*-
"""
@created: Nov 26 2020
@modified: Mar 24 2021
@author: Yoann Pradat

Tests for _load_cln.py module.
"""

import os
import pandas as pd
from pyprism.util import setwd_to_data
from pyprism.data import load_cln

# User functions =======================================================================================================

def test_load_cln_met500():
    df_cln = load_cln(study="met500")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="met500", mode="all")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="met500", multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="met500", identifiers=["MO_1031", "MO_1107"], identifiers_name="Subject_Id")
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2

    df_cln = load_cln(study="met500", identifiers=["MO_1031", "MO_1107"], identifiers_name="Subject_Id",
                      multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2


def test_load_cln_prism():
    df_cln = load_cln(study="prism")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="prism", mode="all")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="prism", multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="prism", identifiers=["201410666GH", "200905404EN"], identifiers_name="Subject_Id")
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2

    df_cln = load_cln(study="prism", identifiers=["201410666GH", "200905404EN"], identifiers_name="Subject_Id",
                      multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2


def test_load_cln_tcga():
    df_cln = load_cln(study="tcga")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="tcga", mode="brca")
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="tcga", multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame

    df_cln = load_cln(study="tcga", identifiers=["TCGA-OR-A5J1", "TCGA-OR-A5J2"], identifiers_name="Subject_Id")
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2

    df_cln = load_cln(study="tcga",identifiers=["TCGA-OR-A5J1", "TCGA-OR-A5J2"], identifiers_name="Subject_Id",
                     multi_cols=True)
    assert type(df_cln) == pd.core.frame.DataFrame
    assert df_cln.shape[0] == 2
