# -*- coding: utf-8 -*-
"""
Created on Aug 18 2021

@author: Yoann Pradat

Tests for _load_summary.py module.
"""

from pyprism.data import load_summary_rna_fus, load_summary_rna_gex, load_summary_wes_mut

def test_load_summary_rna_gex_met500():
    df_sum = load_summary_rna_gex(study="met500", level="genes", metric="counts", test_mode=False)
    df_sum = load_summary_rna_gex(study="met500", level="genes", metric="counts", test_mode=True)
    df_sum = load_summary_rna_gex(study="met500", level="genes", metric="TPM", test_mode=False)
    df_sum = load_summary_rna_gex(study="met500", level="genes", metric="TPM", test_mode=True)


def test_load_summary_rna_fus_met500():
    df_sum = load_summary_rna_fus(study="met500", mode="arriba")
    df_sum = load_summary_rna_fus(study="met500", mode="ericscript")
    df_sum = load_summary_rna_fus(study="met500", mode="pizzly")
    df_sum = load_summary_rna_fus(study="met500", mode="starfusion")


def test_load_summary_rna_gex_prism():
    df_sum = load_summary_rna_gex(study="prism", level="genes", metric="counts", test_mode=False)
    df_sum = load_summary_rna_gex(study="prism", level="genes", metric="counts", test_mode=True)
    df_sum = load_summary_rna_gex(study="prism", level="genes", metric="TPM", test_mode=False)
    df_sum = load_summary_rna_gex(study="prism", level="genes", metric="TPM", test_mode=True)
    df_sum = load_summary_rna_gex(study="prism", level="transcripts", metric="counts", test_mode=False)
    df_sum = load_summary_rna_gex(study="prism", level="transcripts", metric="TPM", test_mode=False)


def test_load_summary_rna_fus_prism():
    df_sum = load_summary_rna_fus(study="prism", mode="arriba")
    df_sum = load_summary_rna_fus(study="prism", mode="ericscript")
    df_sum = load_summary_rna_fus(study="prism", mode="pizzly")
    df_sum = load_summary_rna_fus(study="prism", mode="starfusion")


def test_load_summary_rna_gex_tcga():
    df_sum = load_summary_rna_gex(study="tcga", level="genes", metric="counts", test_mode=False)
    df_sum = load_summary_rna_gex(study="tcga", level="genes", metric="counts", test_mode=True)
    df_sum = load_summary_rna_gex(study="tcga", level="genes", metric="TPM", test_mode=False)
    df_sum = load_summary_rna_gex(study="tcga", level="genes", metric="TPM", test_mode=True)


def test_load_summary_rna_fus_tcga():
    df_sum = load_summary_rna_fus(study="tcga", mode="aggregated")


def test_load_summary_rna_gex_tcga():
    df_sum = load_summary_rna_gex(study="tcga_6_samples", level="genes", metric="counts", other_mode="prism")
    df_sum = load_summary_rna_gex(study="tcga_6_samples", level="genes", metric="counts", other_mode="zheng")
