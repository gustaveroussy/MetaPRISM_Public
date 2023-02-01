# -*- coding: utf-8 -*-
"""
@created: Feb 10 2021
@modified: Apr 25 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions
    load_dsg()
"""

from ._filepaths import _get_filepath_dsg_curated
from ._util_load import load_from_data
from pyprism.util import setwd_to_data

# PRISM ================================================================================================================

def load_dsg(study: str):
    """
    Load the design of the specified study.

    Parameters
    ----------
    study: str
        Only 'prism' is supported

    Returns
    -------
    df_des: DataFrame
    """
    if study.lower() == "prism" or study.lower() == "met500":
        df = load_from_data(_get_filepath_dsg_curated(study))
        if study=="prism":
            # Keep 1 FASTQ per sample
            df = df.dropna(subset=["Datafile_Path"])
            if "Comment_Sample(Y.Pradat)" in df:
                df = df[df["Comment_Sample(Y.Pradat)"].isnull()]
            assert df.shape[0]==df["Sample_Id"].nunique()
        return df
    else:
        raise ValueError("Unsupported value '%s' of study. Choose 'prism' or 'met500'.")
