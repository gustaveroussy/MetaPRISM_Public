# -*- coding: utf-8 -*-
"""
@created: 07 Jul 21
@modified: 16 Feb 22
@authors: Antoine Laine, Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Compare the matrix of signatures scores (ssGSEA) of TCGA samples computed on 2 different versions of the expression
table of TCGA
    1. TCGA expression data used in the paper PMID: 31808800. Expression table available at
    https://stanfordmedicine.app.box.com/s/lu703xuaulfz02vgd2lunxnvt4mfvo3q.
    2. TCGA expression data used in the paper PMID: 34019806. Signatures table available at
    https://github.com/BostonGene/MFP/tree/905fa3cfba58cc379d2af5720235430186863ce6/Cohorts/Pan_TCGA and copied into
    resources/bagaev_2021/Pan_TCGA/signatures.tsv
"""

import argparse
import numpy as np
import pandas as pd
from pyprism.data import load_cln
from scipy.stats import pearsonr, spearmanr, kendalltau

# plot
import seaborn as sns
import matplotlib.pyplot as plt
from prettypy.heatmap import plot_heatmap, HeatmapConfig


# functions ============================================================================================================

def load_bagaev_signatures(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df = df.rename(columns={"Unnamed: 0": "Signature"})
    df = df.set_index(keys="Signature")
    return df


def load_prism_signatures(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df = df.rename(columns={"Unnamed: 0": "Signature"})
    df = df.set_index(keys="Signature")
    return df


def load_bagaev_annotations(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df = df.rename(columns={"Unnamed: 0": "Subject_Id"})
    return df


def load_prism_annotations(filepath):
    df = pd.read_csv(filepath, sep="\t")
    return df


def select_common_cols(df_a, df_b):
    cols_a = df_a.columns.tolist()
    cols_b = df_b.columns.tolist()
    cols_c = list(set(cols_a).intersection(set(cols_b)))
    return df_a[cols_c], df_b[cols_c]


def rename_prism_samples(df_score):
    df_score.columns = list(map(lambda x: x[:12], df_score.columns))
    return df_score


def align_rows_and_cols(df_a_score, df_b_score):
    return df_a_score.loc[df_b_score.index, df_b_score.columns], df_b_score


def pairwise_correlations(df_a, df_b, method):
    names = df_a.columns.tolist()
    corrs = []
    for name in names:
        if method=="pearson":
            corrs.append(pearsonr(df_a[name], df_b[name])[0])
        elif method=="spearman":
            corrs.append(spearmanr(df_a[name], df_b[name])[0])
        elif method=="kendall":
            corrs.append(kendalltau(df_a[name], df_b[name])[0])
        else:
            raise ValueError("Unsupported value of --method.")

    return pd.Series(corrs, index=names)


def normalize_signatures(df):
    return ((df.T - df.T.median())/(df.T.mad())).T


def draw_heatmap(df_correl, method, output):
    config = HeatmapConfig()
    config.figure["figsize"] = (8,8)
    config.figure["n_grid"] = 20
    config.heatmap["ticks_labelsize"] = 6
    config.heatmap["xticks_labelrotation"] = 90
    config.heatmap["yticklabels"] = True
    config.cbar["ticks_labelsize"] = 8
    config.cbar["title_fontsize"] = 10
    config.cbar["title_pad"] = 10
    config.cbar["title"] = method.title()
    config.cbar["xy"] = (0,0.33)
    config.cbar["fraction"] = 0.3
    config.cbar["cmap"] = sns.color_palette("Reds", n_colors=3, as_cmap=True)


    correl_min = np.floor(df_correl.min().min()*100)/100
    correl_max = 1
    config.cbar["boundaries"] = [np.round(c, 2) for c in np.linspace(correl_min, correl_max, num=10)]

    fig, axes = plot_heatmap(df_correl, config)
    title = "%s correlation signatures scores TCGA\nBagaev_2021 vs META-PRISM" % method.title()
    plt.title(title, fontsize=14, loc="left", fontweight="bold", pad=20)
    plt.savefig(output, bbox_inches="tight", dpi=300)


def main(args):
    # samples are in columns, signatures in rows
    df_a_score = load_bagaev_signatures(args.bagaev_scores)
    df_b_score = load_prism_signatures(args.prism_scores)
    df_a_annot = load_bagaev_annotations(args.bagaev_annots)
    df_b_annot = load_prism_annotations(args.prism_annots)

    df_b_score = rename_prism_samples(df_b_score)
    df_a_score, df_b_score = select_common_cols(df_a_score, df_b_score)
    df_a_score, df_b_score = align_rows_and_cols(df_a_score, df_b_score)
    df_a_annot = df_a_annot.loc[df_a_annot["Subject_Id"].isin(df_a_score.columns)]
    df_b_annot = df_b_annot.loc[df_b_annot["Subject_Id"].isin(df_b_score.columns)]

    # compute correlation between vectors of samples for each signature
    s_correls = []
    names = []
    for tcga_project in df_a_annot["TCGA_project"].unique():
        tcga_project_samples = df_a_annot.loc[df_a_annot["TCGA_project"]==tcga_project, "Subject_Id"]
        df_a = df_a_score[tcga_project_samples]
        df_b = df_b_score[tcga_project_samples]
        s_correl = pairwise_correlations(df_a.T, df_b.T, args.method)
        s_correls.append(s_correl)
        names.append(tcga_project)

    df_correl = pd.concat(s_correls, axis=1)
    df_correl.columns = names

    # save
    draw_heatmap(df_correl, args.method, args.output)

# run ==============================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare signature scores.')
    parser.add_argument('--bagaev_scores', type=str, default="resources/bagaev_2021/Pan_TCGA/signatures.tsv",
                        help='Path to Bagaev table of TCGA scores.')
    parser.add_argument('--bagaev_annots', type=str, default="resources/bagaev_2021/Pan_TCGA/annotation.tsv",
                        help='Path to Bagaev table of TCGA annotations.')
    parser.add_argument('--prism_scores', type=str, help='Path to prism table of TCGA annotations.',
                        default="../../../results/immuno_analysis/tcga/tables/signatures_mfp_model.tsv")
    parser.add_argument('--prism_annots', type=str, help="Path to prism table of TCGA annotations.",
                        default="../../../data/tcga/clinical/curated_other/cln_tcga_all_curated.tsv")
    parser.add_argument('--method',type=str, default="pearson",
                        help='Default : pearson. Correlation method to use ("pearson","kendall","spearman").')
    parser.add_argument('--output',  type=str,
                        default="../../../results/immuno_analysis/tcga/checks/mfp_tcga_check_signatures.pdf",
                        help='Path to output graph.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n")

    main(args)
