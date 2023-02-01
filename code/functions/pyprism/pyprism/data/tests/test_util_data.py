# -*- coding: utf-8 -*-
"""
@created: Nov 8 2020
@modified: Nov 8 2020
@author: Yoann Pradat

Tests for _util_data.py module.
"""

import pandas as pd
from pyprism.data import load_wes_mut, preprocess_wes_mut

# User functions =======================================================================================================

def test_preprocess_wes_mut():
    cohort = "prism"
    df_mut = load_wes_mut(study=cohort, mode="oncotator_filtered")

    df_mut_ok = preprocess_wes_mut(df_mut, cohort, cols_bio=["Sample_Id"])
    assert "Sample_Id" in df_mut_ok

    df_mut_ok = preprocess_wes_mut(df_mut, cohort, cols_cln=["Subject_Id", "Project_TCGA_More"])
    assert "Sample_Id" not in df_mut_ok
    assert "Subject_Id" in df_mut_ok
    assert "Project_TCGA_More" in df_mut_ok

    df_mut_ok = preprocess_wes_mut(df_mut, cohort, cols_cln=["Project_TCGA_More"])
    assert "Sample_Id" not in df_mut_ok
    assert "Subject_Id" not in df_mut_ok
    assert "Project_TCGA_More" in df_mut_ok

    df_mut_ok = preprocess_wes_mut(df_mut, cohort, cols_bio=["Biopsy_Date"], cols_cln=["Project_TCGA_More"])
    assert "Sample_Id" not in df_mut_ok
    assert "Biopsy_Date" in df_mut_ok
    assert "Subject_Id" not in df_mut_ok
    assert "Project_TCGA_More" in df_mut_ok

    df_mut_ok = preprocess_wes_mut(df_mut, cohort, cols_cln=["Subject_Id"], select_pairs=True)
    col_tsid_m = "Tumor_Sample_Barcode"
    col_nsid_m = "Matched_Norm_Sample_Barcode"
    col_pid = "Tumor_Vs_Normal"
    df_mut_ok[col_pid] = df_mut_ok[[col_tsid_m, col_nsid_m]].apply("_vs_".join, axis=1)
    assert df_mut_ok[col_pid].nunique()==df_mut_ok["Subject_Id"].nunique()
