# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 2020
@author: Yoann Pradat
Tests for _neighborhood.py module.
"""

from collections import namedtuple
import numpy as np
import pandas as pd
import os
from pyprism.stats import find_intersections_lists
from pyprism.stats import NeighborhoodAnalysis

import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

dir_file = os.path.dirname(os.path.realpath(__file__))

def test_find_intersections_lists():
    # sizes 2,2: one intersection
    l_1 = [(0,1), (1,0)]
    l_2 = [(0,0.75), (1,0.25)]
    assert find_intersections_lists(l_1, l_2) == [(0,0)]

    # sizes 2,2: no intersection
    l_1 = [(0,1), (1,0)]
    l_2 = [(0,0.75), (1,-0.25)]
    assert find_intersections_lists(l_1, l_2) == []

    # sizes >2,2: one intersection
    l_1 = [(0,1), (0.5, 0.5), (1,0)]
    l_2 = [(0,0.75), (1,0.35)]
    assert find_intersections_lists(l_1, l_2) == [(0,0)]

    # sizes >2,2: one intersection
    l_1 = [(0,1), (0.5, 0.5), (1,0)]
    l_2 = [(0,0.75), (1,0.15)]
    assert find_intersections_lists(l_1, l_2) == [(1,0)]

    # sizes >2,2: no intersection
    l_1 = [(0,1), (0.5, 0.5), (1,0)]
    l_2 = [(0,0.75), (1,-0.25)]
    assert find_intersections_lists(l_1, l_2) == []

    # sizes >2,>2: two intersection
    l_1 = [(0,1), (0.5, 0.5), (1,0)]
    l_2 = [(0,0.75), (0.75, 0.65), (1.15, -1)]
    assert find_intersections_lists(l_1, l_2) == [(0,0), (1,1)]

    # sizes >2,>2: two intersection
    l_1 = [(0,1), (0.5, 0.5), (1,0), (3, -3)]
    l_2 = [(0,0.75), (0.75, 0.65), (1.15, -1), (2, -2.5), (4, -2.75)]
    assert find_intersections_lists(l_1, l_2) == [(0,0), (1,1), (2,3)]


def test_neighborhood_analysis_simu():
    Simulation = namedtuple("Simulation", "X y S E m_0 m_1")

    def simulate_signal(F: int, N: int) -> Simulation:
        """
        Simulate a matrix of observations of F variables in N experiments/individuals
        P(X | y = 0) = m_0' * S + E
        P(X | y = 1) = m_1' * S + E
        with
            * m_0, m_1 binary vectors of size (F,1)
            * S a vector of randomly distributed vars of size (F,1) (signal)
            * E a vector of randomly distributed vars of size (F,1) (noise)
        """

        np.random.seed(123)
        y = pd.Series(np.random.randint(low=0, high=2, size=N))

        # masks for the variables that carry signal
        m_0 = np.random.binomial(n=1, p=0.05, size=F)
        m_1 = np.random.binomial(n=1, p=0.05, size=F)

        # drop variables that carry signals for both classes
        m_c = (m_0 > 0) & (m_1 > 0)
        m_c = m_c.astype(int)
        m_0 = m_0 - m_c
        m_1 = m_1 - m_c

        print("%d/%d var with signal to class 1" % (sum(m_1 > 0), F))
        print("%d/%d var with signal to class 0" % (sum(m_0 > 0), F))

        M_0 = np.repeat(m_0.reshape(-1,1), N, axis=1)
        M_1 = np.repeat(m_1.reshape(-1,1), N, axis=1)

        S = pd.DataFrame(np.random.exponential(scale = 1.5, size=(F, N)))
        S_0 = M_0 * S * (1-y)
        S_1 = M_1 * S * y
        E = pd.DataFrame(np.random.randn(F, N))
        X = S_0 + S_1 + E

        return Simulation(X=X, y=y, S=S, E=E, m_0=m_0, m_1=m_1)

    # def view_signal(simu: Simulation, ax_corr, ax_cbar) -> None:
    #     X = simu.X
    #     y = simu.y
    #     m_0 = simu.m_0
    #     m_1 = simu.m_1

    #     rows_order = X[m_1==1].index.tolist() + X[m_0==1].index.tolist()
    #     cols_order = X.loc[:,y==1].columns.tolist() + X.loc[:,y==0].columns.tolist()
    #     X_viz = X.iloc[rows_order, cols_order]

    #     lim_d = X.quantile(0.01).min()
    #     lim_u = X.quantile(0.99).max()

    #     n_colors = 13
    #     lims = np.linspace(lim_d, lim_u, n_colors)
    #     cmap = sns.diverging_palette(240, 10, n=n_colors, l=50, sep=1, as_cmap=True)

    #     sns.heatmap(
    #         X_viz,
    #         linecolor="royalblue",
    #         linewidths=0.75,
    #         cmap=cmap,
    #         vmin=lim_d,
    #         vmax=lim_u,
    #         square=True,
    #         xticklabels=X_viz.columns,
    #         yticklabels=X_viz.index,
    #         ax=ax_corr,
    #         norm=cm.colors.BoundaryNorm(boundaries=lims, ncolors=256),
    #         cbar_ax=ax_cbar,
    #         cbar_kws={
    #             'spacing': 'uniform',
    #             'ticks': lims[::2],
    #             'format': '%.1g',
    #         }
    #     )

    #     ax_corr.tick_params(axis="x", which="both", labelsize=6)
    #     ax_corr.tick_params(axis="y", which="both", labelsize=6)
    #     ax_corr.set_xlabel("")
    #     ax_corr.set_ylabel("")

    #     ax_cbar.tick_params(labelsize=8)
    #     ax_cbar.yaxis.set_label_position("left")
    #     ax_cbar.set_ylabel(ylabel='values', size=10)

    # simulate noisy observations with signal
    F, N = 1000, 20
    simu = simulate_signal(F, N)

    # # visualize signal
    # fig = plt.figure(figsize=(4,9))

    # main_gs  = gridspec.GridSpec(nrows = 10, ncols = 20, wspace = 0, hspace = 0)
    # left_gs  = main_gs[:, :18]
    # right_gs = main_gs[4:7, 18]
    # ax_corr = plt.subplot(left_gs)
    # ax_cbar = plt.subplot(right_gs)

    # view_signal(simu, ax_corr, ax_cbar)
    # plt.savefig(os.path.join(dir_file, "test_neighborhood_simulated_signal.pdf"), bbox_inches="tight")

    def recall_precision(vars_pred, vars_true):
        tp = len(set(vars_pred).intersection(set(vars_true)))
        fp = len(set(vars_pred).difference(set(vars_true)))
        fn = len(set(vars_true).difference(set(vars_pred)))

        recall = tp/(tp+fn)
        precision = tp/(tp+fp)
        return recall, precision

    #### CLASS 1
    # run neighborhood analysis
    neigh = NeighborhoodAnalysis(corr_type="binary", grid_lims=(-0.5, 1.5), grid_size=250, n_perm=5, seed=123, n_jobs=1)
    neigh.fit(simu.X,simu.y)

    # # visualize cumulative curves
    # fig, ax = plt.subplots(figsize=(8,8))
    # neigh.plot(ax, [0.01, 0.05, 0.5], colors=["red", "orange", "gold"])
    # plt.savefig(os.path.join(dir_file, "test_neighborhood_cumulative_curves_class_1.pdf"), bbox_inches="tight")

    # get significant variables to class 1
    vars_pred = neigh.get_significant_variables(alpha=0.01)
    vars_true = simu.X.iloc[simu.m_1 == 1].index.tolist()

    recall, precision = recall_precision(vars_pred, vars_true)
    print("correlation to class 1: recalle %.3g; precision %.3g" % (recall, precision))

    #### CLASS 0
    # run neighborhood analysis
    neigh = NeighborhoodAnalysis(corr_type="binary", grid_lims=(-0.5, 1.5), grid_size=250, n_perm=5, seed=123, n_jobs=1)
    neigh.fit(simu.X,1-simu.y)

    # # visualize cumulative curves
    # fig, ax = plt.subplots(figsize=(8,8))
    # neigh.plot(ax, [0.01, 0.05, 0.5], colors=["red", "orange", "gold"])
    # plt.savefig(os.path.join(dir_file, "test_neighborhood_cumulative_curves_class_0.pdf"), bbox_inches="tight")

    # get significant variables to class 0
    vars_pred = neigh.get_significant_variables(alpha=0.01)
    vars_true = simu.X.iloc[simu.m_0 == 1].index.tolist()

    recall, precision = recall_precision(vars_pred, vars_true)
    print("correlation to class 1: recalle %.3g; precision %.3g" % (recall, precision))

# def test_neighborhood_golub():
#     def get_df_train():
#         folder   = "/Users/ypradat/Documents/data/Broad_Institute/Leukemia/Golub_1999/data"
#         filename = "ALL_vs_AML_train_set_38_sorted.res.txt"
#         filepath = os.path.join(folder, filename)
# 
#         df_train = pd.read_csv(
#             filepath_or_buffer = filepath,
#             sep = "\t"
#         )
# 
#         #### drop null lines and columns
#         df_train = df_train.iloc[2:, :]
#         df_train = df_train.loc[:, df_train.isnull().mean(axis=0) < 1]
# 
#         #### rename columns containin A, P, M entries for each sample
#         old_columns = df_train.columns
#         new_columns = []
# 
#         for i, c in enumerate(old_columns):
#             if c.startswith("Unnamed"):
#                 new_columns.append("type_%s" % old_columns[i-1])
#             else:
#                 new_columns.append(c)
# 
#         df_train.columns = new_columns
# 
#         return df_train
# 
# 
#     def get_X_train():
#         df_train = get_df_train()
# 
#         X_train = df_train[["Accession"] + [x for x in df_train.columns if x.startswith("AML") or x.startswith("ALL")]]
#         X_train = X_train.set_index("Accession")
# 
#         return X_train
# 
# 
#     X_train = get_X_train()
# 
#     #### deal with negative values by offset
#     #### rationale: some beads are less fluorescent than the background controls on the slide
#     #### intensity values are in fact relative to background levels and negative means no detectable
#     #### cDNA for that probe
# 
#     X_train = X_train.clip(lower=1)
#     labels = pd.Series([0 if x.startswith("ALL") else 1 for x in X_train.columns], index=X_train.columns)
# 
#     neigh = NeighborhoodAnalysis(corr_type="binary", grid_lims=(-0.5, 1.75), grid_size=250, n_perm=400, seed=123, n_jobs=6)
#     neigh.fit(X_train,labels)
# 
#     N_P = neigh._N_P
#     C_X = neigh._C_X
# 
# 
#     fig, ax = plt.subplots(figsize=(8,8))
#     for quantile, ls, lw, color in zip([0.01, 0.05, 0.5], ['-','-','-'], [1, 2, 2], ['red', 'red', 'gold']):
#         cdf_q = pd.DataFrame({"u": neigh._grid, "n": N_P.quantile(1-quantile, axis=1)})
#         cdf_q = cdf_q.loc[cdf_q.n > 0]
#         cdf_q = cdf_q.drop_duplicates(subset=["n"], keep="first") # or last ?
# 
#         ax.plot(
#             cdf_q.u,
#             cdf_q.n,
#             ls    = ls,
#             lw    = lw,
#             color = color,
#             label = "%.2g%%" % (100*quantile)
#         )
# 
#     cdf_x = neigh._get_cdf(C_X)
#     ax.scatter(cdf_x.u, cdf_x.n, color="black", s=10)
# 
#     ax.set_xlim(max(neigh._grid), min(neigh._grid))
#     ax.set_ylim(bottom=1)
#     ax.set_yscale("log")
# 
#     ax.legend(loc="lower right", frameon=False, fontsize=14)
#     ax.spines["right"].set_visible(False)
#     ax.spines["top"].set_visible(False)
# 
#     plt.show()
