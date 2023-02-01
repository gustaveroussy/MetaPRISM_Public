# -*- coding: utf-8 -*-
"""
@created: 09/08/21
@modified: 09/08/21
@authors: Antoine Laine, Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Compare expression between PolyA and HybridSelection samples from the same subjects.
"""

import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import sys
import time as time
from tqdm import tqdm
import random

import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from pyprism.data import load_bio, load_rna, load_table, load_summary_rna

Axes = mpl.axes.Axes
Series = pd.core.series.Series
DataFrame = pd.core.frame.DataFrame
Array = np.ndarray

# functions ============================================================================================================

def load_expression(samples_id):
    df_rna = load_rna(study="met500", metric="TPM", identifiers = samples_id,identifiers_name="Sample_Id")
    df_sum = load_summary_rna(study="met500", metric="TPM", test_mode=False)
    new_col_names = {r["Col_Name"]: r["Sample_Id"] for _, r in df_sum.iterrows()}
    df_rna = df_rna.rename(columns=new_col_names)

    return df_rna


def connect_gene_ids_to_symbols(df_rna, df_genes):
    df_rna = df_rna.rename(columns={"ensembl_gene_id": "Ensembl_Gene_Id_Version"})
    df_genes = df_genes[["Ensembl_Gene_Id_Version", "Hugo_Symbol"]]
    df_rna = df_rna.merge(df_genes, how="left", on="Ensembl_Gene_Id_Version")

    return df_rna


def clear_gene_ids(df_rna):
    for col in ["Ensembl_Gene_Id_Version", "Ensembl_Gene_Id"]:
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
    return df_rna


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


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Compare signature scores.')
    parser.add_argument('--gene_table_met500', type=str,
                        default="workflow/results/inputs/gene_tables/gencode_v27_coding_no-mt_no-histones_updated.tsv",
                        help="Path to gene table linking gene ids or old symbols to updated gene symbols.")
    parser.add_argument('--output', type=str,
                        default="workflow/results/plots/met500",
                        help='Path to output folder')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    # run ==============================================================================================================

    ######## PolyA VS Hybrid ########
    df_bio = load_bio(study="met500", identifiers_name="Sample_Type", identifiers=["RNA_T"])
    df_bio_sel = df_bio.loc[df_bio["Library_Selection"].isin(["PolyA", "Hybrid Selection"])]
    df_pat_sel = df_bio_sel.groupby("Subject_Id")[["Library_Selection", "Biopsy_Site"]].nunique().reset_index()
    df_pat_sel = df_pat_sel.loc[(df_pat_sel["Library_Selection"]==2) & (df_pat_sel["Biopsy_Site"]==1)]
    ran_sel_6_subject = random.sample(list(df_pat_sel.Subject_Id),6)
    df_bio_polyA_samples = []
    df_bio_hybrid_samples = []
    for sub_id in ran_sel_6_subject:
        df_bio_polyA_samples.append(list(df_bio[df_bio.Subject_Id==sub_id][df_bio.Library_Selection=="PolyA"].iloc[[0]].Sample_Id)[0])
        df_bio_hybrid_samples.append(list(df_bio[df_bio.Subject_Id==sub_id][df_bio.Library_Selection=="Hybrid Selection"].iloc[[0]].Sample_Id)[0])
    df_bio_met500_samples = list(set(df_bio_polyA_samples + df_bio_hybrid_samples))

    df_rna_met500 = load_expression(df_bio_met500_samples)
    df_genes_met500 = load_table(args.gene_table_met500)
    df_rna_met500= connect_gene_ids_to_symbols(df_rna_met500, df_genes_met500)
    df_rna_met500 = clear_gene_ids(df_rna_met500)
    df_rna_met500 = summarize_duplicated_symbols(df_rna_met500)
    df_rna_met500 = format_expression(df_rna_met500)

    #Shaping both tables to have similar Colnames
    df_rna_polyA = df_rna_met500[df_bio_polyA_samples]
    df_rna_hybrid = df_rna_met500[df_bio_hybrid_samples]
    sample_to_subject = df_bio[df_bio.Sample_Id.isin(df_bio_met500_samples)][["Sample_Id","Subject_Id"]]
    sample_to_subject = sample_to_subject.set_index("Sample_Id")

    df_rna_polyA = df_rna_polyA.rename(columns=sample_to_subject.to_dict()['Subject_Id'])
    df_rna_hybrid = df_rna_hybrid.rename(columns=sample_to_subject.to_dict()['Subject_Id'])

    common_index = df_rna_polyA.index.intersection(df_rna_hybrid.index)
    df_rna_hybrid = df_rna_hybrid.loc[common_index]
    df_rna_polyA = df_rna_polyA.loc[common_index]

    common_columns = df_rna_polyA.columns.intersection(df_rna_hybrid.columns)
    df_rna_hybrid = df_rna_hybrid.loc[:,common_columns]
    df_rna_polyA = df_rna_polyA.loc[:,common_columns]

    plot_samples = df_rna_hybrid.columns.tolist()[:6]
    df_rna_polyA = df_rna_polyA[plot_samples]
    df_rna_hybrid = df_rna_hybrid[plot_samples]

    #Scatter Plot
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9))
    axes = axes.flatten()

    for plot_sample, (i,ax) in zip(plot_samples, enumerate(axes)):
        s_polyA = df_rna_polyA[plot_sample]
        s_hybrid = df_rna_hybrid[plot_sample]
        l_polyA = "PolyA %s" % plot_sample[:12]
        l_hybrid = "Hybrid %s" % plot_sample[:12]
        density_scatter_plot_gex_ab(ax, s_polyA, s_hybrid, l_polyA, l_hybrid)

    plt.subplots_adjust(left=0.1, bottom=0.075, top=0.90)
    plt.suptitle("tpm PolyA vs Hybrid on 6 MET500 samples", fontsize=20, fontweight="medium")
    plt.savefig(os.path.join(args.output,"compare_tpm_polyA_vs_hybrid_met500.png"), bbox_inches="tight", dpi=300)


    ######## Hybrid VS Hybrid ########
    df_bio_sel = df_bio.loc[df_bio["Library_Selection"]=="Hybrid Selection"]
    df_pat_sel = df_bio_sel.groupby("Subject_Id")[["Library_Selection", "Biopsy_Site","Sample_Id"]].nunique().reset_index()
    df_pat_sel = df_pat_sel.loc[(df_pat_sel["Library_Selection"]==1) & (df_pat_sel["Biopsy_Site"]==1) & (df_pat_sel["Sample_Id"]>=2)]
    ran_sel_6_subject = random.sample(list(df_pat_sel.Subject_Id),6)
    df_bio_hybrid1_samples = []
    df_bio_hybrid2_samples = []
    for sub_id in ran_sel_6_subject:
        df_bio_hybrid1_samples.append(list(df_bio[df_bio.Subject_Id==sub_id][df_bio.Library_Selection=="Hybrid Selection"].iloc[[0]].Sample_Id)[0])
        df_bio_hybrid2_samples.append(list(df_bio[df_bio.Subject_Id==sub_id][df_bio.Library_Selection=="Hybrid Selection"].iloc[[1]].Sample_Id)[0])
    df_bio_met500_samples = list(set(df_bio_hybrid1_samples + df_bio_hybrid2_samples))

    df_rna_met500 = load_expression(df_bio_met500_samples)
    df_genes_met500 = load_table(args.gene_table_met500)
    df_rna_met500= connect_gene_ids_to_symbols(df_rna_met500, df_genes_met500)
    df_rna_met500 = clear_gene_ids(df_rna_met500)
    df_rna_met500 = summarize_duplicated_symbols(df_rna_met500)
    df_rna_met500 = format_expression(df_rna_met500)

    #Shaping both tables to have similar Colnames
    df_rna_hybrid1 = df_rna_met500[df_bio_hybrid1_samples]
    df_rna_hybrid2 = df_rna_met500[df_bio_hybrid2_samples]
    sample_to_subject = df_bio[df_bio.Sample_Id.isin(df_bio_met500_samples)][["Sample_Id","Subject_Id"]]
    sample_to_subject = sample_to_subject.set_index("Sample_Id")

    df_rna_hybrid1 = df_rna_hybrid1.rename(columns=sample_to_subject.to_dict()['Subject_Id'])
    df_rna_hybrid2 = df_rna_hybrid2.rename(columns=sample_to_subject.to_dict()['Subject_Id'])

    common_index = df_rna_hybrid1.index.intersection(df_rna_hybrid2.index)
    df_rna_hybrid2 = df_rna_hybrid2.loc[common_index]
    df_rna_hybrid1 = df_rna_hybrid1.loc[common_index]

    common_columns = df_rna_hybrid1.columns.intersection(df_rna_hybrid2.columns)
    df_rna_hybrid2 = df_rna_hybrid2.loc[:,common_columns]
    df_rna_hybrid1 = df_rna_hybrid1.loc[:,common_columns]

    plot_samples = df_rna_hybrid2.columns.tolist()[:6]
    df_rna_hybrid1 = df_rna_hybrid1[plot_samples]
    df_rna_hybrid2 = df_rna_hybrid2[plot_samples]

    #Scatter Plot
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,9))
    axes = axes.flatten()

    for plot_sample, (i,ax) in zip(plot_samples, enumerate(axes)):
        s_hybrid1 = df_rna_hybrid1[plot_sample]
        s_hybrid2 = df_rna_hybrid2[plot_sample]
        l_hybrid1 = "hybrid1 %s" % plot_sample[:12]
        l_hybrid2 = "hybrid2 %s" % plot_sample[:12]
        density_scatter_plot_gex_ab(ax, s_hybrid1, s_hybrid2, l_hybrid1, l_hybrid2)

    plt.subplots_adjust(left=0.1, bottom=0.075, top=0.90)
    plt.suptitle("tpm Hybrid vs Hybrid on 6 MET500 samples", fontsize=20, fontweight="medium")
    plt.savefig(os.path.join(args.output,"compare_tpm_hybrid_vs_hybrid_met500.png"), bbox_inches="tight", dpi=300)




    #Distribution plot
    polyA = df_rna_polyA.iloc[:, 0]
    hybrid = df_rna_hybrid.iloc[:, 0]
    filter = (polyA + hybrid)
    filter_index = filter[filter>5].index
    diff = (polyA - hybrid)/(polyA + hybrid)
    diffnona = diff.dropna()
    diff_final = diffnona.loc[filter_index].sort_values(ascending=True)
    #diff_final.to_csv(os.path.join(args.output,"met500_polyA_vs_hybrid_expression.tsv"),sep="\t")
    sns.displot(diff_final.values)
    plt.savefig(os.path.join(args.output,"PolyAvsHybrid.png"))

    polyA = df_rna_polyA.iloc[:, 0]
    polyA_2 = df_rna_polyA.iloc[:, 1]
    filter = (polyA + polyA_2)
    filter_index = filter[filter>5].index
    diff = (polyA - polyA_2)/(polyA + polyA_2)
    diffnona = diff.dropna()
    diff_final = diffnona.loc[filter_index].sort_values(ascending=True)
    #diff_final.to_csv(os.path.join(args.output,"met500_polyA_vs_polyA_expression.tsv"),sep="\t")
    sns.displot(diff_final.values)
    plt.savefig(os.path.join(args.output,"PolyAvsPolyA.png"))
