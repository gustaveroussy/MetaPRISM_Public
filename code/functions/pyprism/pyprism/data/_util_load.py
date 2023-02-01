# -*- coding: utf-8 -*-
"""
@created: Dec 15 2020
@modified: Apr 27 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Helper functions for loading.

Functions
    def check_filepath_exists_in_data(filepath) -> bool:
    def delineate_data(df, col_name, prefix="", suffix="", format=["Dates", "Cycle", "Dose"], delineate_dates=False
                       convert_dates_to_datetime=True, sep=","):
    def load_table(path: str, sep=None, header=0, dtype=None, **kwargs) -> DataFrame:
    def load_from_data(path: str, sep=None, header=0, dtype=None, **kwargs) -> DataFrame:
    def read_header(path: str) -> list:
    def subset_data(df: DataFrame, values: Iterable[Any]=None, col_name: str=None):
"""

import os
import gzip
import numpy as np
import pandas as pd
from pandas.api.types import is_datetime64_any_dtype
import re
from pyprism.util import setwd_to_data, explode_df
from typing import Iterable, Any, Callable

DataFrame = pd.core.frame.DataFrame


def convert_dates_to_str(df, strftime="%d/%m/%Y"):
    df_c = df.copy()
    for x in df.columns:
        if is_datetime64_any_dtype(df_c[x]):
            df_c.loc[:,x] = df_c[x].apply(lambda x: x.strftime(strftime) if not pd.isnull(x) else "").astype(str)
    return df_c


def check_filepath_exists_in_data(filepath) -> bool:
    if filepath is None:
        filepath_exists = False
    else:
        cwd = setwd_to_data()
        filepath_exists = os.path.exists(filepath)
        os.chdir(cwd)

    return filepath_exists


def delineate_data(df, col_name, prefix="", suffix="", format=["Dates", "Cycle", "Dose"], delineate_dates=False,
                   convert_dates_to_datetime=True, sep=",", name="Name", remove=False, preserve_index=False, **kwargs):
    """
    Delineate treatments from compact format to more user-friendly format.

    Parameters
    ----------
    df: DataFrame
        Dataframe that must contain the column `col_name`.
    col_name: str
        Name of the column to be delineated
    prefix: str, default=""
        Prefix to the new columns that contain the delineated data.
    suffix: str, default=""
        Suffix to the new columns that contain the delineated data.
    format: str, default=["Dates", "Cycle", "Dose"]
        List of delineated column names corresponding the format of the compressed data.
    delineate_dates: str, default=True
        Set to True only if "Dates" is in `format`. Split Dates into Date_Min and Date_Max.
    convert_dates_to_datetime: str, default=True
        Used only if delineate_dates=True. Should dates be converted to datetime type?
    sep: str, default=","
        The delimiter
    name: str, default="Name"
        Name given to the column containing the value before the brackets. Prefix and suffixes specified will be added
        to this name.
    remove: bool, default=False,
        If True, remove the column `col_name`.
    preserve_index: bool, default=False,
        If True, preserve the original index names.
    **kwargs:
        Additional arguments passed to pd.to_datetime function. Used only if convert_dates_to_datetime is True.
    """
    df_c = df.copy()
    df_c[col_name] = df_c[col_name].fillna("Unknown")
    df_c = explode_df(df_c, cols=[col_name], sep="|", preserve_index=preserve_index)

    s_extract = df_c[col_name].str.extract(r"(?<=\()(.*)(?=\))")
    df_c_split = s_extract[0].str.split(sep, expand=True)
    df_c_split.columns = ["%s%s%s" % (prefix, f, suffix) for f in format][:df_c_split.shape[1]]

    if delineate_dates and "Dates" in format:
        col_dates = "%sDates%s" % (prefix,suffix)
        col_date_min = "%sDate_Min%s" % (prefix,suffix)
        col_date_max = "%sDate_Max%s" % (prefix,suffix)
        df_c_split_dates = df_c_split[col_dates].str.split("_to_", expand=True)

        # no value with "_to_" tag
        if df_c_split_dates.shape[1] == 1:
            df_c_split[[col_date_min, col_date_max]] = np.nan
        else:
            df_c_split[[col_date_min, col_date_max]] = df_c_split[col_dates].str.split("_to_", expand=True)

    if convert_dates_to_datetime:
        cols_date = [x for x in df_c_split.columns if "Date" in x and not "Dates" in x]
        for col_date in cols_date:
            df_c_split[col_date] = pd.to_datetime(df_c_split[col_date], **kwargs)

    c_name = "%s%s%s" % (prefix,name,suffix)
    df_c[c_name] = df_c[col_name].apply(lambda x: re.sub(r"\([^)]*\)", "", x))

    df_c[col_name] = df_c[col_name].replace("Unknown", np.nan)
    df_c[c_name] = df_c[c_name].replace("Unknown", np.nan)

    df_c = pd.concat((df_c, df_c_split), axis=1)

    if remove:
        del df_c[col_name]

    return df_c


def _read_excel_with_dates(path, skiprows=None, header=0, dtype=None, cols_date=None, dayfirst=True, format=None, **kwargs):
    if dtype is None and cols_date is not None:
        dtype = {k: "str" for k in cols_date}
    elif dtype is not None and cols_date is not None:
        for col in cols_date:
            dtype[col] = "str"
    elif cols_date is not None:
        raise ValueError("If using both cols_date and dtype arguments, please specify a dictionary for dtype.")

    df_cln = pd.read_excel(path, na_values=["na", "na ", "NA", ".", "not available"], header=header,
                           skiprows=skiprows, parse_dates=False, engine="openpyxl", dtype=dtype, **kwargs)

    if cols_date is not None:
        for col in cols_date:
            if format is not None:
                df_cln.loc[:, col] = pd.to_datetime(df_cln[col], format=format, errors="coerce")
            else:
                df_cln.loc[:, col] = pd.to_datetime(df_cln[col], dayfirst=dayfirst, errors="coerce")

    return df_cln



def load_table(path: str, header=0, dtype=None, header_prefix=None, **kwargs) -> DataFrame:
    if header_prefix is not None:
        skiprows = len(read_header(path=path, prefix=header_prefix))
    else:
        skiprows = None

    if path.endswith("xlsx"):
        df = _read_excel_with_dates(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)
    elif ".txt" in path or ".tsv" in path or ".maf" in path:
        if "sep" not in kwargs:
            df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, sep="\t", **kwargs)
        else:
            df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)
    else:
        df = pd.read_csv(path, skiprows=skiprows, header=header, dtype=dtype, **kwargs)

    return df


def load_from_data(path: str, header=0, dtype=None, header_prefix=None, **kwargs) -> DataFrame:
    cwd = setwd_to_data()
    df = load_table(path, header, dtype, header_prefix, **kwargs)
    os.chdir(cwd)

    return df


def make_column_groups(df: DataFrame, cols2groups) -> DataFrame:
    df_new = df.loc[:, cols2groups.keys()].copy()
    tuples =  [(v,k) for k, v in cols2groups.items()]
    df_new.columns = pd.MultiIndex.from_tuples(tuples, names=["col_group", "col_name"])
    return df_new


def save_df(df, path, sheet_name="Sheet1", dates_to_str=True, index=False, verbose=True, **kwargs):
    if dates_to_str:
        df_save = convert_dates_to_str(df)
    else:
        df_save = df

    if type(path) == pd.io.excel._xlsxwriter.XlsxWriter or path.endswith("xlsx"):
        df_save.to_excel(path, sheet_name=sheet_name, engine="openpyxl", index=index, **kwargs)
    elif path.endswith(".txt") or path.endswith(".tsv"):
        df_save.to_csv(path, index=index, sep="\t", **kwargs)
    else:
        df_save.to_csv(path, index=index, **kwargs)

    if verbose:
        print("-file saved at %s" % path, flush=True)


def read_header(path: str, prefix: str) -> list:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith(prefix)]
    return header


def save_df_to_data(df, path, sheet_name="Sheet1", dates_to_str=True, index=False, verbose=True, **kwargs):
    cwd = setwd_to_data()
    save_df(df, path, sheet_name, dates_to_str, index, verbose, **kwargs)
    os.chdir(cwd)


def set_dtype_columns(df, cols, dtype):
    cols_c = list(set(df.columns).intersection(set(cols)))
    for col_c in cols_c:
        df[col_c] = df[col_c].astype(dtype)
    if dtype=="str":
        df = df.replace("nan", np.nan)
    return df


def subset_data(df: DataFrame, values: Iterable[Any]=None, col_name: str=None):
    if values is None or col_name is None:
        return df
    else:
        return df.loc[df[col_name].isin(values),:]
