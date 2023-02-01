# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 20:29 2019

@author: Yoann Pradat

The present file implements usual statistical tests used in the analysis of tables.
"""

import numpy as np
import pandas as pd

from scipy.stats import chi2, norm, fisher_exact
from typing import Union, List, Tuple

# type alias
Vector = Union[List[float], np.ndarray, pd.core.series.Series]
Table  = Union[List[List[float]], np.ndarray, pd.core.frame.DataFrame]

def cont_table(A: Vector, B: Vector) -> Table:
    """
    Compute 2 x 2 contignecy table from A and B binary vectors.

    Parameters
    ----------
    A: array-like
        A vector of binary entries
    B: array-like
        A vector of binary entries

    Returns
    -------
    tab: array-like
        A 2x2 contigency table.
    """
    tab = np.zeros((2, 2))
    A_anti = np.where(A==1, 0, 1)
    B_anti = np.where(B==1, 0, 1)
    tab[0,0] = np.sum(A*B)
    tab[0,1] = np.sum(A*B_anti)
    tab[1,0] = np.sum(A_anti*B)
    tab[1,1] = np.sum(A_anti*B_anti)
    return tab

def odds_ratio(tab: Table) -> float:
    """
    Computes the odds ratio of a contigency table
    -------------------
        a        b
        c        d
    -------------------
    as (a/b)/(c/d) or ad/bc

    Parameters
    ----------
    tab: array-like
        The table.

    Returns
    -------
    odds_ratio: float
    """

    if tab[0,1] == 0 or tab[1,0] == 0:
        odds_ratio = tab[0,0] * tab[1,1] / max(tab[1,0], tab[0,1], 1)
    else:
        odds_ratio = tab[0,0] * tab[1,1] / (tab[1,0] * tab[0,1])

    return odds_ratio

def fisher_test(tab: Table) -> Tuple[float, float]:
    """
    Fisher exact test

    Parameters
    ----------
    tab: array-like
        The table.

    Returns
    -------
    odds_ratio, p_value: float, float
        The odds ratio and the p-value.
    """
    odds_ratio, p_value = fisher_exact(tab)
    return odds_ratio, p_value


def table_test(tab: Table, type: str='MH') -> Tuple[float, int, float]:
    """
    Function for testing independence hypotheses in tables.

    Mantel-Haenszel statistic
        # In this test we assume that margins of the table are fixed.
        # If the table is:
        #           Event
        # Group     Yes   No
        #     0     d0    n0-d0   n0
        #     1     d1    n1-d1   n1
        #            d     n-d    n
        #
        # The random variable D0 ~ hypergeometric law. From this we deduce
        # expected value and variance. The statistic is then the usual Wald
        # test
        #
        # Null hypothesis H0: No association between event and group

    The Pearson statistic
        # Assess equality of the 2 groups. It is approximately equivalent
        # to the Mantel-Haenszel statistic but assumes only fixed row margins.
        # The statistic is
        # ---           2
        # \        (o-e)
        # /        -----
        # ---        e
        #
        # Where the sum is over all cells -> nrows*ncols terms
        #
        # Null hypothesis H0: No association between event and group
    Parameters
    ----------
    tab: array-like
        The table.
    type: str
        The name of the statistical test desired.

    Returns
    -------
    statistic, df, p_value: float, int, float
        The value of the statistic, the number of degrees of freedom and the associated p-value.
    """
    tab = np.array(tab)

    if type == 'MH' and tab.shape == (2, 2):
        # degrees of freedom
        df = (tab.shape[0] - 1) * (tab.shape[1] - 1)
        n0 = tab.sum(axis=1)[0]
        d0 = tab[0, 0]
        d = tab.sum(axis=0)[0]
        n = tab.sum(axis=1).sum(axis=0)
        # expected value
        E = n0 * d / n
        # variance
        V = n0 * (n - n0) * d * (n - d) / (n ** 2 * (n - 1))

        statistic = (d0 - E) ** 2 / V
        p_value = 1 - chi2.cdf(statistic, df)

    elif type == 'MH' and tab.shape != (2, 2):
        raise ValueError("Mantel-Haenszel statistic is not yet implemented for tables greater than 2x2.")

    elif type == 'Pearson':
        # degrees of freedom
        df = (tab.shape[0] - 1) * (tab.shape[1] - 1)
        n = tab.sum(axis=1).sum(axis=0)
        cols_margins = tab.sum(axis=0)[:, np.newaxis]
        rows_margins = tab.sum(axis=1)[:, np.newaxis]
        # expected table
        tab_exp = np.dot(rows_margins, cols_margins.T) / n

        statistic = ((tab - tab_exp) ** 2 / tab_exp).sum(axis=1).sum(axis=0)
        p_value = 1 - chi2.cdf(statistic, df)

    else:
        raise ValueError("Unsupported type %s. Choose one of: 'MH', 'Pearson'." % type)

    return statistic, df, p_value
