# -*- coding: utf-8 -*-
"""
@created: Dec 15 2020
@modified: Dec 15 2020
@author: Yoann Pradat

Tests for _load_bio.py module.
"""

import pandas as pd
from pyprism.data import load_bio

# User functions =======================================================================================================

def test_load_bio_met500():
    df_bio = load_bio(study="met500")
    assert type(df_bio) == pd.core.frame.DataFrame

    df_bio = load_bio(study="met500", mode="all")
    assert type(df_bio) == pd.core.frame.DataFrame

def test_load_bio_prism():
    df_bio = load_bio(study="prism")
    assert type(df_bio) == pd.core.frame.DataFrame

    df_bio = load_bio(study="prism", mode="all")
    assert type(df_bio) == pd.core.frame.DataFrame

def test_load_bio_tcga():
    df_bio = load_bio(study="tcga")
    assert type(df_bio) == pd.core.frame.DataFrame

    df_bio = load_bio(study="tcga", mode="all")
    assert type(df_bio) == pd.core.frame.DataFrame

    df_bio = load_bio(study="prism", identifiers=["201410666GH", "200905404EN"], identifiers_name="Subject_Id")
    assert type(df_bio) == pd.core.frame.DataFrame

    df_bio = load_bio(study="tcga", identifiers=["TCGA-02-0047", "TCCA-02-0055"], identifiers_name="Subject_Id")
    assert type(df_bio) == pd.core.frame.DataFrame
