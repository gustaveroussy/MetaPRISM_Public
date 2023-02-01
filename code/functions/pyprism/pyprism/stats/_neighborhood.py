# -*- coding: utf-8 -*-
"""
Created on Tue May 12 2020

@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Neighborhood analysis as described in Golub et al. 1999. 10.1126/science.286.5439.531
    _are_intersecting_lines
    _are_intersectable_lists
    _find_intersections_lists
    _is_decreasing_list_points

    find_intersections_lists
    NeighborhoodAnalysis
"""

from   bisect import bisect
import concurrent.futures
from joblib import Parallel, delayed
import os
import time
from typing import Callable, Union

import matplotlib as mlp
import matplotlib.pyplot as plt
import numpy  as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

Series = pd.core.series.Series
DataFrame = pd.core.frame.DataFrame
Array = np.ndarray
Vector = Union[list, Series, Array]
Axes = mlp.axes.Axes

def _are_intersecting_lines(a_1: float, b_1: float, a_2: float, b_2: float) -> bool:
    """
    Assess whether two lines intersect within the limits of the points provided.
    """
    # check that xcoords are in increasing order
    assert a_1[0] < b_1[0]
    assert a_2[0] < b_2[0]

    m_1 = (b_1[1] - a_1[1])/(b_1[0] - a_1[0])
    m_2 = (b_2[1] - a_2[1])/(b_2[0] - a_2[0])

    p_1 = - m_1 * a_1[0] + a_1[1]
    p_2 = - m_2 * a_2[0] + a_2[1]

    if m_1 == m_2:
        # parallel lines
        return False
    else:
        x_i = (p_2 - p_1)/(m_1 - m_2)
        if x_i >= max(a_1[0], a_2[0]) and x_i <= min(b_1[0], b_2[0]):
            return True
        else:
            return False

def _are_intersectable_lists(l_1: list, l_2: list) -> bool:
    """
    Assess if two lists of points in decreasing order (i.e decreasing second coordintate) may intersect.
    """
    if l_1[-1][1] >= l_2[0][1] or l_2[-1][1] >= l_1[0][1]:
        # lists are completely separable by a horizontal line
        return False
    elif l_1[-1][0] <= l_2[0][0] or l_2[-1][0] <= l_1[0][0]:
        # lists are completely separable by a vertical line
        return False
    else:
        return True


def _find_intersections_lists(l_1: list, l_2: list, i_1: int, i_2: int, p: list) -> None:
    """
    Recursively find all intersection points of two lists of points sorted in decreasing order.
    """
    if _are_intersectable_lists(l_1, l_2):
        if len(l_1) == 2 and len(l_2) == 2:
            if _are_intersecting_lines(l_1[0], l_1[-1], l_2[0], l_2[-1]):
                p.append((i_1, i_2))
        else:
            if len(l_2) == 2:
                if len(l_1) % 2 == 1:
                    _find_intersections_lists(l_1[:(len(l_1)+1)//2], l_2, i_1, i_2, p)
                    _find_intersections_lists(l_1[(len(l_1)-1)//2:], l_2, i_1 + (len(l_1)-1)//2, i_2, p)
                else:
                    _find_intersections_lists(l_1[:len(l_1)//2], l_2, i_1, i_2, p)
                    _find_intersections_lists(l_1[len(l_1)//2-1:], l_2, i_1 + len(l_1)//2-1, i_2, p)
            else:
                if len(l_2) % 2 == 1:
                    _find_intersections_lists(l_1, l_2[:(len(l_2)+1)//2], i_1, i_2, p)
                    _find_intersections_lists(l_1, l_2[(len(l_2)-1)//2:], i_1, i_2 + (len(l_2)-1)//2, p)
                else:
                    _find_intersections_lists(l_1, l_2[:len(l_2)//2], i_1, i_2, p)
                    _find_intersections_lists(l_1, l_2[len(l_2)//2-1:], i_1, i_2 + len(l_2)//2-1, p)

def _is_decreasing_list_points(l: list) -> bool:
    l_x = [e[0] for e in l]
    l_y = [e[1] for e in l]
    return (l_x == sorted(l_x)) & (l_y == sorted(l_y)[::-1])

def find_intersections_lists(l_1: list, l_2: list) -> list:
    """
    Find coordinates (in list indices) of the intersection points of two lists of points sorted in decreasing order.
    """
    assert _is_decreasing_list_points(l_1)
    assert _is_decreasing_list_points(l_2)

    points = []
    _find_intersections_lists(l_1, l_2, 0, 0, points)
    return points

class NeighborhoodAnalysis(object):
    """
    Compute cumulative count curves of variables correlation to
        - a vector of observed labels
        - vectors of random permutations of the observed labels

    Using a discretized grid over specified bounds, the number of variables that have a correlation with observed
    labels and the random permutations of the labels higher than a fixed valued from the grid are computed.
    From the distribution of the number of variables above a fixed value in permutations, one may estimate which
    variables have a correlation to the labels that is "higher than chance" where the chance level can be chosen.
    """
    def __init__(self, corr_type: str, corr_func: Callable=None, grid_lims: tuple=(-0.5, 1.5), grid_size: int=100, n_perm: int=500, seed: int=0, n_jobs: int=1, verbose: int=1):
        """
        Parameters
        ----------
        corr_type: str
            A string among the following:
                - "custom": use the argument corr_func
                - "binary": labels should be exactly 0 or 1. Returns the ratio of the difference between the mean of
                    class 1 and the mean of class 0 to the sum of the standard deviations of class 1 and class 0.
                - "spearman": Spearman correlation coefficient.
                - "pearson": Pearson correlation coefficient.
        corr_func: callable
            A callable function that takes exactly two arguments: a series of variable values and a series of
            labels and returns a float value. Used only if corr_type is set to "custom".
        grid_lims: tuple
            Limits of the grid that define the cumulative count curves
        grid_size: int
            Number of unique values where the cumulative counts are computed.
        n_perm: int
            Number of random permutations of the labels.
        seed: int
            Seed of numpy random generator.
        n_jobs: int
            Number of jobs used by joblib Parallel function.
        verbose : int
            Verbosity level of jolib Parallel function.
        """
        if corr_type == "custom":
            self.corr_func = corr_func
        else:
            self.corr_func = lambda v,c: self._corr_func(v,c,corr_type)

        self.grid_lims = grid_lims
        self.grid_size = grid_size
        self.n_perm = n_perm
        self.seed = seed
        self.n_jobs = n_jobs
        self.verbose = verbose
        self._intersection_points = None

    def fit(self, X: DataFrame, y: Vector):
        """
        Parameters
        ----------

        X: DataFrame
            Data of size (F,N) organized in variables x observations.
        y: Vector
            Numpy array, pandas series or list of size (N,) defining the observations labels.
        """
        # matrix of labels permutations
        self._labels_perm = self._get_labels_permutations(y, self.seed)

        # vector of correlations between X rows and labels
        self._C_X = self._correlations_to_labels(X, y)

        # vector of cumulative counts
        self._cdf_x = self._get_cdf(self._C_X)

        # matrix of correlations between X rows and labels permutations
        self._C_P = self._correlations_to_labels_permutations(X, self._labels_perm)

        # matrix of variable counts with correlation above values of a grid
        self._N_P = self._compute_N_P(self._C_P)

        return self

    def _get_labels_permutations(self, labels: Series, seed: int) -> list:
        """
        Compute list of random permutations of the labels.
        """
        np.random.seed(seed)
        labels_perm = [np.random.permutation(labels) for _ in range(self.n_perm)]
        return labels_perm

    def _corr_func(self, v: Vector, c: Vector, corr_type:str) -> float:
        c = np.array(c)
        v = np.array(v)
        if corr_type == "binary":
            assert set(np.unique(c)).issubset(set([0,1]))

            mu_0 = np.mean(v[c==0])
            mu_1 = np.mean(v[c==1])
            sig_0 = np.std(v[c==0])
            sig_1 = np.std(v[c==1])

            if sig_0 + sig_1 == 0:
                corr = 0
            else:
                corr = (mu_1 - mu_0)/(sig_0 + sig_1)
        elif corr_type == "pearson":
            corr, p = pearsonr(v, c)
        elif corr_type == "spearman":
            corr, p = spearmanr(v, c)
        else:
            raise ValueError("Invalid corr_type: %s. Choose one of 'binary', 'pearson' or 'spearman'" % corr_type)

        return corr

    def _correlations_to_labels(self, X: DataFrame, labels: Series) -> Series:
        """
        Compute vector of correlations between each row of X and labels.
        """
        return X.apply(lambda v: self.corr_func(v, labels), axis=1)

    def _correlations_to_labels_permutations(self, X: DataFrame, labels_perm: list) -> DataFrame:
        """
        Compute correlation between each row of the matrix X and each permutation of the vector of labels.
        The correlations are stored in a dataframe of size (F, n_perm).
        """
        parallel = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)
        C_P = parallel(delayed(self._correlations_to_labels)(X, labels) for labels in labels_perm)
        return pd.DataFrame(pd.concat(C_P, axis=1))

    def _compute_N_P(self, C_P: DataFrame) -> DataFrame:
        """
        Compute the number of variables that have a corrleation above certain values of a grid.
        The numbers are stored in a dataframe of size (grid_size, n_perm).
        """
        self._grid = np.linspace(self.grid_lims[0], self.grid_lims[1], self.grid_size)
        N_P = np.zeros((self.grid_size, self.n_perm))

        for t in range(self.grid_size):
            for p in range(self.n_perm):
                N_P[t,p] = (C_P.iloc[:,p] >= self._grid[t]).sum()

        return pd.DataFrame(N_P)

    def _get_cdf(self, C: Series) -> DataFrame:
        us = C.sort_values().unique()
        ns = [(C >= u).sum() for u in us]
        return pd.DataFrame({"u": us, "n": ns})


    def _get_intersection_points_xq(self, quantile: float=0.01) -> Array:
        cdf_q = pd.DataFrame({"u": self._grid, "n": self._N_P.quantile(quantile, axis=1)})
        cdf_q = cdf_q.loc[cdf_q.n > 0]
        cdf_q = cdf_q.drop_duplicates(subset=["n"], keep="first") # or last ?
        cdf_q = cdf_q.reset_index(drop=True)

        l_x = list(self._cdf_x.to_records(index=False))
        l_q = list(cdf_q.to_records(index=False))

        points = find_intersections_lists(l_x, l_q)

        return points, cdf_q

    def get_significant_variables(self, alpha: float=0.01) -> Array:
        """
        For a given alpha, a "chance curve" is built by computing the 1-alpha quantile of the distribution of counts
        at each value of the grid. Significant variables are the the variables that are in the part of the observed
        curve that is above the chance curve.

        Parameters
        ----------
        alpha: float
            Level of "significance" i.e quantile that is will be used to decided upon the "chance curve" to be used to
            compare against the observed curve.
        """
        points, cdf_q = self._get_intersection_points_xq(quantile = 1 - alpha)

        # get the abscisse just after the last intersection points
        u_lim  = self._cdf_x.iloc[points[-1][0]+1]["u"]
        n_lim  = self._cdf_x.iloc[points[-1][0]+1]["n"]

        vars_at_alpha = self._C_X.loc[self._C_X >= u_lim].index.tolist()
        assert len(vars_at_alpha) == n_lim

        return self._C_X.loc[vars_at_alpha]

    def plot(self, ax: Axes, alphas: list=[0.01, 0.05, 0.5], lss: list=["-", "-", "-"], lws: list=[2,2,2], colors: list=["dodgerblue", "red", "orange"], hlines: list=[True, True, True], legend_fontsize: int=10) -> None:

        ax.scatter(self._cdf_x.u, self._cdf_x.n, color="black", s=10)

        for alpha, ls, lw, color, hline in zip(alphas, lss, lws, colors, hlines):
            points, cdf_q = self._get_intersection_points_xq(quantile = 1 - alpha)

            ax.plot(
                cdf_q.u,
                cdf_q.n,
                ls    = ls,
                lw    = lw,
                color = color,
                label = "%.2g%%" % (100*alpha)
            )

            if hline:
                n_lim  = self._cdf_x.iloc[points[-1][0]+1]["n"]
                u_lim  = self._cdf_x.iloc[points[-1][0]+1]["u"]
                ax.hlines(
                    y     = n_lim,
                    xmin  = u_lim,
                    xmax  = max(self._grid),
                    ls    = "--",
                    lw    = 1,
                    color = color
                )
                ax.annotate(text="%d" % n_lim, xy=(max(self._grid), n_lim), xycoords="data", color=color)

        ax.set_xlim(max(self._grid), min(self._grid))
        ax.set_ylim(bottom=1)
        ax.set_yscale("log")

        ax.legend(loc="lower right", frameon=False, fontsize=legend_fontsize)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
