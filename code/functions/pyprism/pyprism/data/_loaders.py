# -*- coding: utf-8 -*-
"""
@modified: Dec 22 2020
@created: Jun 19 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Code for the class
    LoaderCln
"""

from itertools import chain
from abc import ABC

import numpy as np
import pandas as pd

from pyprism.data import load_bio, load_cln, load_summary_rna_gex

DataFrame = pd.core.frame.DataFrame()

class LoaderUtils(ABC):
    def _drop_nas(self, df_cln: DataFrame, cln_names: list, verbose: bool=True) -> DataFrame:
        keep = df_cln[cln_names].isnull().sum(axis=1) == 0
        if verbose:
            mess = "subset %d/%d samples with complete clnotations:" % (sum(keep), len(keep))
            inde = " " * len(mess) + " "
            clns = ("\n%s" % inde).join(cln_names)
            print(mess, clns, flush=True)

        return df_cln.loc[keep,:]

    def _cln_keep(self, df_cln: DataFrame, cln2keep: dict, verbose: bool=True) -> DataFrame:
        keep_tot = pd.Series([True]*df_cln.shape[0], index=df_cln.index)
        if cln2keep is not None:
            for cln_name, cln_values in cln2keep.items():
                if type(cln_values) != list:
                    cln_values = [cln_values]
                keep = df_cln.loc[keep_tot, cln_name].isin(cln_values)
                if verbose:
                    mess = "subset %d/%d samples on the following value(s) of %s:" % (sum(keep), sum(keep_tot), cln_name)
                    inde = " " * len(mess) + " "
                    clns = ("\n%s" % inde).join([str(val) for val in cln_values])
                    print(mess, clns, flush=True)
                keep_tot = keep & keep_tot

        return df_cln.loc[keep_tot,:]


    def _cln_apply(self, df_cln: DataFrame, cln2apply: dict, verbose: bool=True) -> DataFrame:
        keep_tot = pd.Series([True]*df_cln.shape[0], index=df_cln.index)
        if cln2apply is not None:
            for cln_cols, (cln_name, cln_apply) in cln2apply.items():
                keep = df_cln[[x for x in cln_cols]].apply(cln_apply, axis=1)
                if verbose:
                    mess = "subset %d/%d samples on the following filter: %s" % (sum(keep), sum(keep_tot), cln_name)
                    print(mess, flush=True)
                keep_tot = keep & keep_tot

        return df_cln.loc[keep_tot,:]


    def _cln_rename(self, df_cln: DataFrame, old2new_cln: dict) -> DataFrame:
        df_cln = df_cln.rename(columns=old2new_cln)
        return df_cln

    def _init_dict_with_keys(self, dt: dict, keys: list) -> dict:
        if dt is None:
            return {key: None for key in keys}
        else:
            return {key: dt[key] if key in dt.keys() else None for key in keys}


class LoaderCln(LoaderUtils):
    """
    Class for loading clnotations (i.e clinical) data with specific settings.
    """
    def __init__(self, studies: list, cln2keep_all: dict=None, cln2keep_each: dict=None,
                 cln2apply_each: dict=None, cln2apply_all: dict=None,  verbose: bool=True):
        """
        Parameters
        ----------
        studies: list
            Studies names
        cln2keep_each: dictionary
            Values to subset the data on for study-specific clnotations (e.g 'Primary_Doubt' in 'prism')
        cln2keep_all: dictionary
            Values to subset the data on for study-shared clnotations (e.g Age)
        cln2apply_each: dict
            Callable functions used subset the data on for study-specific clnotations (e.g 'Primary_Doubt' in 'prism')
        cln2apply_all: dictionary
            Callable functions used to subset the data on for study-shared clnotations (e.g Age)
        verbose: bool, optional
            For intermediate messages.

        Returns
        -------
        None

        """
        self.studies = studies
        self.cln2keep_each = self._init_dict_with_keys(cln2keep_each, studies)
        self.cln2keep_all = cln2keep_all
        self.cln2apply_each = self._init_dict_with_keys(cln2apply_each, studies)
        self.cln2apply_all = cln2apply_all
        self.verbose = verbose

    def load(self) -> DataFrame:
        studies = self.studies
        cln2keep_each = self.cln2keep_each
        cln2keep_all = self.cln2keep_all
        cln2apply_each = self.cln2apply_each
        cln2apply_all = self.cln2apply_all

        dfs_cln = []
        for study in self.studies:
            if study == "prism":
                df_cln = self._make_cln_prism(cln2keep=self.cln2keep_each[study],
                                              cln2apply=self.cln2apply_each[study])
            elif study == "tcga":
                df_cln = self._make_cln_tcga(cln2keep=self.cln2keep_each[study],
                                             cln2apply=self.cln2apply_each[study])
            df_cln.insert(0, "Study", study)
            dfs_cln.append(df_cln)

        cols_all = [set(df_cln.columns) for df_cln in dfs_cln]
        cols_com = list(set.intersection(*cols_all))
        df_cln = pd.concat(dfs_cln, axis=0).reset_index(drop=True)

        if self.cln2keep_all is not None:
            df_cln = self._drop_nas(df_cln, list(self.cln2keep_all.keys()), self.verbose)
            df_cln = self._cln_keep(df_cln, self.cln2keep_all, self.verbose)

        if self.cln2apply_all is not None:
            df_cln = self._cln_apply(df_cln, self.cln2apply_all, self.verbose)

        return df_cln

    def _make_cln_prism(self, cln2keep: dict=None, cln2apply: dict=None) -> DataFrame:
        if self.verbose:
            print("aggregating clnotations for PRISM ...")

        old2new_cln = {"Age_At_Biopsy": "Age"}

        df_cln = load_cln(study="prism", multi_cols=False)
        df_cln = self._cln_rename(df_cln, old2new_cln)

        if cln2keep is not None:
            df_cln = self._drop_nas(df_cln, cln2keep.keys(), self.verbose)
            df_cln = self._cln_keep(df_cln, cln2keep, self.verbose)

        if cln2apply is not None:
            df_cln = self._drop_nas(df_cln, list(set(chain(*cln2apply.keys()))), self.verbose)
            df_cln = self._cln_apply(df_cln, cln2apply, self.verbose)

        return df_cln


    def _make_cln_tcga(self, cln2keep: dict=None) -> DataFrame:
        if self.verbose:
            print("aggregating clnotations for TCGA ...")

        old2new_cln = {"Age_At_Diagnosis": "Age",
                       "Tumor_Sample": "Patient"}

        df_cln = load_cln(study="tcga", multi_cols=False)
        df_cln = self._cln_rename(df_cln, old2new_cln)

        if cln2keep is not None:
            df_cln = self._drop_nas(df_cln, cln2keep.keys(), self.verbose)
            df_cln = self._cln_keep(df_cln, cln2keep, self.verbose)

        if cln2apply is not None:
            df_cln = self._drop_nas(df_cln, list(set(chain(*cln2apply.keys()))), self.verbose)
            df_cln = self._cln_apply(df_cln, cln2apply, self.verbose)

        return df_cln
