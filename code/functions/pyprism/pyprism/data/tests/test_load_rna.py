# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 2020

@author: Yoann Pradat

Tests for _load_rna.py module.
"""

from pyprism.data import load_rna_fus, load_rna_gex

def test_load_rna_gex_met500():
    df_rna = load_rna_gex(study="met500", metric="counts", test_mode=True)
    del df_rna

    df_rna = load_rna_gex(study="met500", metric="TPM", test_mode=True, identifiers=["MO_1031", "MO_1107"],
                      identifiers_name="Subject_Id")
    del df_rna


def test_load_rna_fus_met500():
    df_rna = load_rna_fus(study="met500", mode="arriba")
    df_rna = load_rna_fus(study="met500", mode="pizzly")
    df_rna = load_rna_fus(study="met500", mode="oncokb_civic")


def test_load_rna_gex_prism():
    df_rna = load_rna_gex(study="prism", metric="TPM", test_mode=True)
    del df_rna

    df_rna = load_rna_gex(study="prism", metric="TPM", test_mode=True, identifiers=["201410666GH", "200905404EN"],
                      identifiers_name="Subject_Id")
    del df_rna


def test_load_rna_fus_prism():
    df_rna = load_rna_fus(study="prism", mode="ericscript")
    df_rna = load_rna_fus(study="prism", mode="starfusion")
    df_rna = load_rna_fus(study="prism", mode="oncokb_civic")


def test_load_rna_gex_tcga():
    df_rna = load_rna_gex(study="tcga", metric="TPM", test_mode=True)
    del df_rna

    df_rna = load_rna_gex(study="tcga", metric="TPM", test_mode=True,
                      identifiers=["TCGA-02-0047-01A-01R-1849-01", "TCGA-3X-AAV9-01A-72R-A41I-07"],
                      identifiers_name="Sample_Id")
    del df_rna


def test_load_rna_fus_tcga():
    df_rna = load_rna_fus(study="tcga", mode="aggregated")


def test_load_rna_gex_tcga_6_samples():
    df_rna = load_rna_gex(study="tcga_6_samples", metric="TPM", other_mode="prism")
    del df_rna

    df_rna = load_rna_gex(study="tcga_6_samples", metric="counts", other_mode="zheng")
    del df_rna
