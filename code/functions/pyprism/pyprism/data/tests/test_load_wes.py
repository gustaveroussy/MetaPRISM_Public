# -*- coding: utf-8 -*-
"""
@created: Jan 24 2021
@modified: Aug 30 2021
@author: Yoann Pradat

Tests for _load_wes_mut.py module.
"""

import os
import pandas as pd
from pyprism.data import load_wes_mut

def test_load_wes_mut_met500():
    df_maf = load_wes_mut(study="met500", mode="somatic_maf")
    df_maf = load_wes_mut(study="met500", mode="somatic_filters")
    df_maf = load_wes_mut(study="met500", mode="somatic_oncokb")
    df_maf = load_wes_mut(study="met500", mode="somatic_civic")

def test_load_wes_mut_prism():
    df_maf = load_wes_mut(study="prism", mode="somatic_maf")
    df_maf = load_wes_mut(study="prism", mode="somatic_filters")
    df_maf = load_wes_mut(study="prism", mode="somatic_oncokb")
    df_maf = load_wes_mut(study="prism", mode="somatic_civic")

def test_load_wes_mut_tcga():
    # df_maf = load_wes_mut(study="tcga", mode="somatic_maf")
    # df_maf = load_wes_mut(study="tcga", mode="somatic_filters")
    df_maf = load_wes_mut(study="tcga", mode="somatic_oncokb")
    df_maf = load_wes_mut(study="tcga", mode="somatic_civic")
