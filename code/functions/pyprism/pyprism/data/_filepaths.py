# -*- coding: utf-8 -*-
"""
@created: 19/04/21
@modified: 16/06/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Define filepaths to raw and curated data tables.
"""

from pyprism.util import setwd_to_data

import os
import numpy as np
import pandas as pd
import yaml

def _get_filepaths():
    cwd = setwd_to_data()
    filepath = "organisation/DATA_ORGANISATION.yaml"
    FILEPATHS = yaml.load(open(filepath, "r"), Loader=yaml.FullLoader)
    os.chdir(cwd)
    return FILEPATHS

def _get_filepath_bio_curated(study, mode):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["bio"][mode]


def _get_filepath_cln_curated(study, mode):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["cln"][mode]


def _get_filepath_dsg_curated(study):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["dsg"]


def _get_filepath_ids_curated(study):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["ids"]


def _get_filepath_res(database, name):
    FILEPATHS = _get_filepaths()
    return FILEPATHS["resources"][database][name]


def _get_filepath_colors() -> str:
    FILEPATHS = _get_filepaths()
    return "data_overview/colors/colors.xlsx"


def _get_filepath_wes_mut(study, mode):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["wes_mut"][mode]


def _get_filepath_rna_fus(study, mode):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["rna_fus"][mode]


def _get_filepath_rna_gex(study, level="genes", metric="counts", test_mode=False, other_mode=None):
    FILEPATHS = _get_filepaths()
    if level not in FILEPATHS[study]["rna_gex"]:
        if other_mode is not None:
            mode = other_mode
        else:
            if test_mode:
                mode = "test"
            else:
                mode = "full"
        return FILEPATHS[study]["rna_gex"][mode][level][metric]["data"]
    else:
        return FILEPATHS[study]["rna_gex"][level][metric]["data"]


def _get_filepath_summary_rna_fus(study, mode="unfiltered"):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["rna_fus"]["summary"][mode]


def _get_filepath_summary_rna_gex(study, level="genes", metric="counts", test_mode=False, other_mode=None):
    FILEPATHS = _get_filepaths()
    if level not in FILEPATHS[study]["rna_gex"]:
        if other_mode is not None:
            mode = other_mode
        else:
            if test_mode:
                mode = "test"
            else:
                mode = "full"
        return FILEPATHS[study]["rna_gex"][mode][level][metric]["summary"]
    else:
        return FILEPATHS[study]["rna_gex"][level][metric]["summary"]


def _get_filepath_summary_wes_mut(study, mode):
    FILEPATHS = _get_filepaths()
    return FILEPATHS[study]["wes_mut"]["summary"][mode]
