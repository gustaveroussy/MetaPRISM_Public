# -*- coding: utf-8 -*-
"""
@created: 29/07/21
@modified: 29/07/21
@authors: Antoine Laine, Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Compare expression of Xena to prism for TCGA samples.
"""

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import time as time
from tqdm import tqdm

import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from pyprism.data import load_table

Axes = mpl.axes.Axes
Series = pd.core.series.Series
DataFrame = pd.core.frame.DataFrame
Array = np.ndarray

# functions ============================================================================================================


def load_expression(expression_data, expression_summary=None, **kwargs):
    df_rna = pd.read_table(expression_data, **kwargs)

    if expression_summary is not None:
        df_sum = pd.read_table(expression_summary)
        df_sum["Match_Id"] = df_sum["Sample_Id"].str[:15]
        df_sum = df_sum.drop_duplicates(subset=["Match_Id"], keep=False)

        cols_keep = ["ensembl_gene_id", "geneName", "geneID"] + df_sum["Col_Name"].tolist()
        df_rna = df_rna[[x for x in cols_keep if x in df_rna]]
        old2new = {r["Col_Name"]: r["Match_Id"] for _, r in df_sum.iterrows()}
        df_rna = df_rna.rename(columns=old2new)

    return df_rna


def connect_gene_ids_to_symbols(df_rna, df_genes):
    df_rna = df_rna.rename(columns={"ensembl_gene_id": "Ensembl_Gene_Id_Version", "geneID": "Ensembl_Gene_Id_Version"})
    df_genes = df_genes[["Ensembl_Gene_Id_Version", "Hugo_Symbol"]]
    df_rna = df_rna.merge(df_genes, how="left", on="Ensembl_Gene_Id_Version")

    return df_rna


def clear_gene_ids(df_rna):
    for col in ["Ensembl_Gene_Id_Version", "Ensembl_Gene_Id", "geneName", "geneID"]:
        if col in df_rna.columns:
            del df_rna[col]

    return df_rna


def summarize_duplicated_symbols(df_rna):
    if df_rna["Hugo_Symbol"].nunique() < df_rna.shape[0]:
        return df_rna.groupby("Hugo_Symbol").mean().reset_index()
    else:
        return df_rna


def format_expression(df_rna):
    df_rna = df_rna.set_index("Hugo_Symbol")
    df_rna.columns = list(map(lambda x: x[:15], df_rna.columns))
    return df_rna


def preprocess_expression(df_rna):
    return np.log2(df_rna+1)


def depreprocess_function(x):
    y = 2**x - 1e-3
    return np.where(np.abs(y) < 1e-5, 0, y)


def depreprocess_expression(df_rna):
    return df_rna.apply(depreprocess_function)


def regression_xy(x: Array, y: Array) -> tuple:
    rho = np.corrcoef(x, y)[0,1]
    X = np.concatenate((np.ones((x.shape[0],1)),x.reshape(-1,1)), axis=1)

    beta = np.linalg.inv((X.T.dot(X))).dot(X.T).dot(y)
    r2 = 1 - np.sum(np.square((y - X.dot(beta))))/np.sum(np.square((y - y.mean())))

    return rho, beta, r2


def density_scatter_plot_gex_ab(ax: Axes, s_a: Series, s_b: Series, l_a: str, l_b:str, lims: tuple=None,
                                log_scale:bool=True) -> None:
    x = s_a.values
    y = s_b.values

    # density colors: long for large datasets
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    z = (z - z.min())/(z.max()-z.min())
    cmap = cm.get_cmap("viridis")
    c = cmap(z)

    # scatter plot
    ax.scatter(x, y, c=c, s=20)
    ax.set_xlabel(l_a, fontsize=15, fontweight="medium")
    ax.set_ylabel(l_b, fontsize=15, fontweight="medium")

    # regression line
    rho, beta, r2 = regression_xy(x,y)
    x_line = np.array([min(x), max(x)])
    y_line = beta[0] + beta[1] * x_line

    label = "LS regression\n  y = %.3g x + %.3g \nCorrelation\n  %.3g\nR-squared\n  %.3g" % (beta[1], beta[0], rho, r2)
    ax.annotate(label, xy=(0.1, 0.6), xycoords='axes fraction')

    # esthetics
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(axis="both", which="both", labelsize=12)

    lims = (1, 1.05*max(max(x), max(y)))
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")


def scatter_plot_gex_ab(ax: Axes, s_a: Series, s_b: Series, l_a: str, l_b:str, log_scale: bool=True) -> None:
    x = s_a.values
    y = s_b.values

    ax.scatter(x, y, s=50, edgecolor='')
    ax.set_xlabel(l_a, fontsize=15, fontweight="medium")
    ax.set_ylabel(l_b, fontsize=15, fontweight="medium")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(axis="both", which="both", labelsize=12)

    lims = (1, 1.05*max(max(x), max(y)))
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    if log_scale:
        ax.set_xscale("log")
        ax.set_yscale("log")



if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Compare signature scores.')
    parser.add_argument('--expression_prism', type=str, help="Path to prism TCGA expression table.")
    parser.add_argument('--summary_prism', type=str, help="Path to prism TCGA summary expression table.")
    parser.add_argument('--gene_table_prism', type=str,
                        default="../../../results/immuno_analysis/gene_tables/gencode_v27_updated.tsv",
                        help="Path to gene table linking gene ids or old symbols to updated gene symbols.")
    parser.add_argument('--expression_xena', type=str,
                        default="resources/ucsc_xena/Pan_TCGA/tcga_target_no_normal_rsem_gene_tpm.gz",
                        help="""Path to Xena TCGA expression table.""")
    parser.add_argument('--gene_table_xena', type=str,
                        default="../../../results/immuno_analysis/gene_tables/gencode_v23_updated.tsv",
                        help="Path to gene table linking gene ids or old symbols to updated gene symbols.")
    parser.add_argument('--output', type=str,
                        default="workflow/results/plots/tcga/compare_tpm_xena_vs_prism.png",
                        help='Path to output table of signature scores.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    # run ==============================================================================================================
    df_rna_xena = load_expression(args.expression_xena, usecols=np.arange(0,12), sep="\t")
    df_genes_xena = load_table(args.gene_table_xena)
    df_rna_xena = df_rna_xena.rename(columns={"sample": "ensembl_gene_id"})
    df_rna_xena = connect_gene_ids_to_symbols(df_rna_xena, df_genes_xena)
    df_rna_xena = clear_gene_ids(df_rna_xena)
    df_rna_xena = summarize_duplicated_symbols(df_rna_xena)
    df_rna_xena = format_expression(df_rna_xena)
    df_rna_xena = depreprocess_expression(df_rna_xena)

    df_rna_prism = load_expression(args.expression_prism, args.summary_prism)
    df_genes_prism = load_table(args.gene_table_prism)
    df_rna_prism = connect_gene_ids_to_symbols(df_rna_prism, df_genes_prism)
    df_rna_prism = clear_gene_ids(df_rna_prism)
    df_rna_prism = summarize_duplicated_symbols(df_rna_prism)
    df_rna_prism = format_expression(df_rna_prism)
    # df_rna_prism = preprocess_expression(df_rna_prism)

    common_index = df_rna_xena.index.intersection(df_rna_prism.index)
    df_rna_prism = df_rna_prism.loc[common_index]
    df_rna_xena = df_rna_xena.loc[common_index]

    common_columns = df_rna_xena.columns.intersection(df_rna_prism.columns)
    df_rna_prism = df_rna_prism.loc[:,common_columns]
    df_rna_xena = df_rna_xena.loc[:,common_columns]

    plot_samples = df_rna_prism.columns.tolist()[:6]
    df_rna_xena = df_rna_xena[plot_samples]
    df_rna_prism = df_rna_prism[plot_samples]

    # plot
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9))
    axes = axes.flatten()

    for plot_sample, (i,ax) in zip(plot_samples, enumerate(axes)):
        s_xena = df_rna_xena[plot_sample]
        s_prism = df_rna_prism[plot_sample]
        l_xena = "Xena %s" % plot_sample[:12]
        l_prism = "Prism %s" % plot_sample[:12]
        density_scatter_plot_gex_ab(ax, s_xena, s_prism, l_xena, l_prism)

    plt.subplots_adjust(left=0.1, bottom=0.075, top=0.90)
    plt.suptitle("tpm Xena vs PRISM on 6 TCGA samples", fontsize=20, fontweight="medium")
    plt.savefig(args.output, bbox_inches="tight", dpi=300)
