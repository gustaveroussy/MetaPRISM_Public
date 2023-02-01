# -*- coding: utf-8 -*-
"""
@created: 27/07/21
@modified: 27/07/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Plot 2d representation of high-dimensional data with flexible coloring.
"""

import argparse
import numpy as np
import pandas as pd

import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

Axes = mpl.axes.Axes
Series = pd.core.series.Series
DataFrame = pd.core.frame.DataFrame
Array = np.ndarray

def run_pca(X: DataFrame, n_components=100, return_CTR=False, return_C02=False):
    """
    Paramters
    ---------
        X: array
            Centered array in format n_samples x n_features
        n_components: int, default=100
            Number of principal components
        CTR: bool, default=False
            If True, computes the contribution of each sample to each component.
        C02: bool, default=False
            If True, computes the squared cosinus (quality of the representation) of each sample on each component.
    """
    n_samples = X.shape[0]
    n_features = X.shape[1]
    pca = PCA(n_components=n_components)

    # center
    pca = pca.fit(X)
    C = pca.transform(X)
    C = pd.DataFrame(C, index=X.index)

    def _get_table_ctr(C: DataFrame) -> DataFrame:
        """
        Table of size n_samples x n_components that gives the contributions of each sample to the construction of
        each of the principal axes. Contribution to u_k principal axis by x_n is:
                (x_n' u_k)^2 / (N * lambda_k)

        x_n'u_k is given by C[n,k] and lambda_k by eigenval[k]

        CTR[n,k] represents the contribution of sample n to the construction of principal axis u_k.
        """
        N = C.shape[0]
        l = (C ** 2).sum(axis=0).values / N
        CTR = (C ** 2) / (N*l)
        return CTR

    def _get_table_co2(C, X):
        """
        Table of size n_samples x n_components that gives the quality of the representation of each sample on each of
        the principal axes. Quality is given by
                (x_n' u_k)^2 / ||y_n||^2

        x_n'u_k is given by C[n,k]

        CO2[n,k] represents the quality of the representation of sample n on the principal axis u_k.
        """
        n = (X ** 2).sum(axis=1)
        C02 = ((C ** 2).T / n).T
        return C02

    if return_CTR and return_C02:
        CTR = _get_table_ctr(C)
        C02 = _get_table_co2(C, X)
        return C, CTR, C02
    elif return_CTR:
        CTR = _get_table_ctr(C)
        return C, CTR
    elif return_C02:
        C02 = _get_table_co2(C, X)
        return C, C02
    else:
        return C


def run_tsne(X: DataFrame):
    """
    Paramters
    ---------
        X: array
            Array in format n_samples x n_features
    """

    tsne = TSNE(
        n_components = 2,
        perplexity   = 30,
        learning_rate = 200,
        n_iter = 1000,
        metric = "euclidean",
        verbose = 1,
        method = "barnes_hut",
        angle = 0.5
    )

    C_pca = run_pca(X, n_components=min(50, X.shape[1]))
    C_tsne = tsne.fit_transform(C_pca)
    C_tsne = pd.DataFrame(C_tsne, index=C_pca.index)

    return C_tsne


class Representation2D(object):
    def __init__(self, C, Y, filepath, title="", prefix_axis_name="Principal", c_x=0, c_y=1,
                 colors=None, markers=None, s=50, figsize=(14,10)):
        """
        Parameters
        ----------
        C: array
            An array of format n_samples x n_components.
        Y: array
            An array of format n_samples x 2 where the 2 columns are Label_Color and Label_Marker
        filepath: str
            Path to where the plot shall be saved.
        title: str
            Title of the plot.
        prefix_axis_name: str
            Prefix to the axes names.
        c_x: int
            Index of the component in x-axis.
        c_y: int
            Index of the component in y-axis.
        colors: dict
            Keys should match values in Y.Label_Color. Values should be colors.
        markers: dict
            Keys should match values in Y.Label_Marker. Values should be colors. May be NULL.
        s: int
            Size of the dots.
        figsize: tuple
            Dimensions of the plot in inches.
        """
        self.C = C
        self.Y = Y
        self.filepath = filepath
        self.title = title
        self.prefix_axis_name = prefix_axis_name
        self.c_x = c_x
        self.c_y = c_y
        self.colors = colors
        self.markers = markers
        self.s = s
        self.figsize = figsize

    def _plot_scatter(self, ax, C, y, colors, alpha: float=1, marker: str='o', s: int=50, edgecolors: str='none'):
        for label in y.unique():
            ax.scatter(
                x          = C.loc[y==label].iloc[:,0],
                y          = C.loc[y==label].iloc[:,1],
                color      = colors[label],
                alpha      = alpha,
                marker     = marker,
                lw         = 0.2,
                s          = s,
                edgecolors = edgecolors
            )

    def _draw_legend(self, ax: Axes, colors: dict=None, markers: dict=None, fontsize: int=12):
        """
        Draw manual legend using two dictionaries.
        """
        handles = []
        if colors is not None:
            for label, color in sorted(colors.items(), key=lambda item: item[0]):
                handle = Line2D([0], [0], marker='s', color=color, label=label, markersize=10, linestyle='None')
                handles.append(handle)
        if markers is not None:
            for label, marker in markers.items():
                handle = Line2D([0], [0], marker=marker, color='black', label=label, markersize=10,
                                markerfacecolor='None', linestyle='None')
                handles.append(handle)

        ax.legend(handles=handles, loc="best", frameon=False, fontsize=fontsize)


    def _draw_spines(self, ax: Axes, label_x: str, label_y: str, fontsize: int):
        """
        Draw beautiful spines.
        """
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_position("zero")
        ax.spines["left"].set_position("zero")

        # x-axis
        ax.set_xlabel("")
        ax.annotate(label_x, size=fontsize, xy=(1, 0), xycoords=('axes fraction', 'data'), xytext=(0, 0),
                    textcoords='offset points', ha='right', va='bottom')
        # y-axis
        ax.set_ylabel("")
        ax.annotate(label_y, size=fontsize, xy=(0, 1), xycoords=('data', 'axes fraction'), xytext=(0, 0),
                    textcoords='offset points', ha='right', va='top', rotation='vertical')

        ax.set_xticklabels([])
        ax.set_yticklabels([])


    def _save_plot(self, filepath: str):
        plt.savefig(filepath, bbox_inches="tight")
        plt.close()
        print("-plot saved at %s" % filepath)


    def draw(self):
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=self.figsize)
        ax.set_title(self.title, fontsize=20)

        if self.markers is None:
            self._plot_scatter(
                ax         = ax,
                C          = self.C.iloc[:,[self.c_x, self.c_y]],
                y          = self.Y["Label_Color"],
                colors     = self.colors,
                marker     = 's',
                s          = self.s,
                edgecolors = 'black'
            )
        else:
            for label, marker in self.markers.items():
                mask_marker = self.Y["Label_Marker"]==label
                self._plot_scatter(
                    ax         = ax,
                    C          = self.C.loc[mask_marker].iloc[:,[self.c_x, self.c_y]],
                    y          = self.Y.loc[mask_marker, "Label_Color"],
                    colors     = self.colors,
                    marker     = marker,
                    s          = self.s,
                    edgecolors = 'black'
                )

        self._draw_legend(ax, colors=self.colors, markers=self.markers, fontsize=12)
        self._draw_spines(ax, label_x="%s Axis %d" % (self.prefix_axis_name, self.c_x + 1),
                          label_y="%s Axis %d" % (self.prefix_axis_name, self.c_y + 1), fontsize=16)
        self._save_plot(self.filepath)
