# -*- coding: utf-8 -*-
"""
@created: Jul 08 2021
@modified: Aug 06 2021

@author: Yoann Pradat

Tests for _load_res.py module.
"""

import os
from pyprism.data import load_resource

def test_load_resource():
    # CIViC
    df_gene = load_resource(database="civic", name="gene")

    # CGC
    df_gene = load_resource(database="cosmic", name="gene")

    # curated
    df_gene = load_resource(database="curated", name="gene")

    # Gencode
    df_gene = load_resource(database="gencode", name="gencode_v23")

    # OncoKB
    df_gene = load_resource(database="oncokb", name="gene")
