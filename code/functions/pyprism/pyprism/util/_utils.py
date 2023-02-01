# -*- coding: utf-8 -*-
"""
@created: 01/24/21
@modified: 01/24/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Useful functions.
"""

import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

def explode_df(df: DataFrame, cols: list, sep=',', fill_value: str='', preserve_index: bool=False) -> DataFrame:
    """
    Expand dataframe entries of the columns specified in l_cols and for which there are multiple values.

    Parameters
    ---------
    df: DataFrame
        Input dataframe on which expansion is performed
    cols: list
        List of columns where expansion is required
    sep: char
        Character separating the multiple values
        Default : ','
    fill_value: bool
        Entry in exanpded dataframe for empty lists
        Default : ''
    preserve_index: bool
        Whether original index should be preserved or not. If set to True, the index of the expanded DataFrame
        will be redundant.
        Default : False

    Returns
    -------
    df_expanded: DataFrame
        Returns  dataframe where entries of any of the columns in l_cols with multiple values have been expanded.
    """
    # transform comma-separated to list
    df = df.assign(**{col:df[col].str.split(sep) for col in cols}).copy()
    if (cols is not None and len(cols) > 0 and not isinstance(cols, (list, tuple, np.ndarray, pd.Series))):
        cols = [cols]
    # calculate lengths of lists
    lens = df[cols[0]].str.len()
    # format NaN to [NaN] and strip unwanted characters
    for col in cols:
        df.loc[df[col].isnull(), col] = df.loc[df[col].isnull(), col].apply(lambda x: [np.nan])
        df.loc[lens > 1, col] = df.loc[lens > 1, col].apply(lambda x: [y.strip() for y in x])
    # all columns except `cols`
    idx_cols = df.columns.difference(cols)
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    df_expanded = (pd.DataFrame({col:np.repeat(df[col].values, lens) for col in idx_cols},
                index=idx).assign(**{col:np.concatenate(df.loc[lens>0, col].values) for col in cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        df_expanded = (df_expanded.append(df.loc[lens==0, idx_cols], sort=False).fillna(fill_value))
    # revert the original index order
    df_expanded = df_expanded.sort_index()
    # reset index if requested
    if not preserve_index:
        df_expanded = df_expanded.reset_index(drop=True)
    return df_expanded


def rm_duplicates_list(seq):
    """
    Remove duplicates in a list while preserving order.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
