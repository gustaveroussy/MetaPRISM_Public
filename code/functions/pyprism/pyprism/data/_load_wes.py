# -*- coding: utf-8 -*-
"""
@created: Jan 24 2021
@modified: Apr 26 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions
    def _build_load_stat_for_wes_mut(df_maf, patient_field, tumor_field, normal_field):
    def load_wes_mut(study: str, identifiers: Iterable[str]=None, identifiers_name: str=None, mode:
        str="somatic_filtered"):
"""

import os
import pandas as pd
from typing import Iterable

from pyprism.util import setwd_to_data
from ._util_load import subset_data, load_from_data
from ._filepaths import _get_filepath_wes_mut

DataFrame = pd.core.frame.DataFrame

def _build_load_stat_for_wes_mut(df_maf, patient_field, tumor_field=None, normal_field=None):
    if tumor_field is None and normal_field is not None:
        df_stat = df_maf.loc[:, [patient_field, normal_field]].drop_duplicates()
        df_stat = df_stat.reset_index(drop=True)
        df_stat_gby = df_stat.groupby(patient_field)
        df_stat = pd.DataFrame({"N_Loaded": df_stat_gby[normal_field].apply(lambda x: x.nunique()),
                                normal_field: df_stat_gby[normal_field].apply(lambda x: "|".join(x))})
        df_stat = df_stat.reset_index(drop=False)

        message = "normal sample"
    else:
        df_stat = df_maf.loc[:, [patient_field, tumor_field, normal_field]].drop_duplicates()
        df_stat = df_stat.reset_index(drop=True)
        tumor_vs_normal = "%s_vs_%s" % (tumor_field, normal_field)
        df_stat.loc[:, tumor_vs_normal] = df_stat[[tumor_field, normal_field]].agg(lambda x: "_vs_".join(x), axis=1)
        df_stat = df_stat[[patient_field, tumor_vs_normal]].reset_index(drop=True)
        df_stat_gby = df_stat.groupby(patient_field)
        df_stat = pd.DataFrame({"N_Loaded": df_stat_gby[tumor_vs_normal].apply(lambda x: x.nunique()),
                                tumor_vs_normal: df_stat_gby[tumor_vs_normal].apply(lambda x: "|".join(x))})
        df_stat = df_stat.reset_index(drop=False)

        message = "couple(s) tumor vs normal"

    # print
    s_n_loaded = df_stat["N_Loaded"].value_counts()

    for n in sorted(s_n_loaded.index):
        if n <= 2:
            print("%d %s loaded for %d patient(s)" % (n, message, s_n_loaded[n]))
        else:
            print("%d %s loaded for %d patient(s)" % (n, message, s_n_loaded[n]))
            print("\t" + "\n\t".join(list(df_stat.loc[df_stat.N_Loaded == n, patient_field].values)))

    return df_stat

# User function ========================================================================================================

def load_wes_mut(study: str, identifiers: Iterable[str]=None, identifiers_name: str=None, mode: str="somatic_maf",
             return_load_stat: bool=False):
    """
    Load mutation tables for the specified study.

    Parameters
    ----------
    study: str
        Choose 'met500', 'tcga' or 'prism'.
    identifiers: Iterable[str], optional.
        If not None, return annotations only for the identifiers values.
    identifiers_name: Iterable[str], optional.
        Name of the column of identifiers.
    mode: str, default="somatic_maf"
        1. If "somatic_maf", load data/[study]/wes/somatic_maf/somatic_calls.maf.gz
        2. If "somatic_filters", load data/[study]/wes/somatic_maf/somatic_calls_filters.tsv.gz
        3. If "somatic_oncokb", load data/[study]/wes/somatic_maf/somatic_calls_oncokb.maf.gz
        4. If "somatic_civic", load data/[study]/wes/somatic_maf/somatic_calls_civic.maf.gz

    return_load_stat: bool, default=True
        If True, returns a tuple of dataframe with first the maf table and second a table providing details about what
        was loaded. If False, returns only the MAF table.


    Returns
    -------
    df_maf, df_stat: tuple of DataFrame
    """

    # load maf file
    print("-reading aggregated MAF file ... ", end="")
    df_maf = load_from_data(_get_filepath_wes_mut(study, mode), header_prefix="##",
                            sep="\t", dtype={"Chromosome": str}, low_memory=False)
    print("ok!")

    # # use metadata to add sample and subject identifiers
    # df_meta = load_summary_wes_mut(study, mode, use_cache=True)
    # del df_meta["N_Mutations"]
    # cols_com = list(set(df_meta.columns).intersection(set(df_maf.columns)))
    # df_maf = df_maf.merge(df_meta, how="left", on=cols_com)

    # subset if any
    df_maf = subset_data(df_maf, values=identifiers, col_name=identifiers_name)

    if return_load_stat:
        if "oncotator" in mode or "mc3" in mode:
            # build a table with one line per subject and the number of pairs tumor-normal loaded for that patient
            df_stat = _build_load_stat_for_wes_mut(df_maf, "Subject_Id",
                                               tumor_field="Tumor_Sample_Id", normal_field="Normal_Sample_Id")
        elif "annovar" in mode:
            df_stat = _build_load_stat_for_wes_mut(df_maf, "Subject_Id", normal_field="Sample_Id")

        return df_maf, df_stat
    else:
        return df_maf
