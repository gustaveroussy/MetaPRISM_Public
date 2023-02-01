# -*- coding: utf-8 -*-
"""
@created: 05 Apr 2022
@modified: 13 Apr 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Double heatmap showing the coocurrence of treatment resistances.
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys
import seaborn as sns

# pyprism
from pyprism.data import load_table, load_cln, load_bio, split_targeted_therapy_targets, load_colors
from pyprism.util import explode_df

# prettypy
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from prettypy.heatmap import build_double_heatmap, plot_double_heatmap, DoubleHeatmapConfig

def get_table_of_drug_data(cohort, tumor_type, level, selection={}):
    col_tt = "Project_TCGA_More"

    # load clinical and biospecimen data
    df_cln = load_cln(cohort)
    df_bio = load_bio(cohort)

    # add columns from biospecimen data
    df_cln["Biopsy_Id"] = df_cln["Biopsy_Selected"].apply(lambda x: x.split("|")[0])
    cols_bio = ["Biopsy_Id", "Biopsy_Site", "Biopsy_Subsite", "Study"]
    df_cln = df_cln.merge(df_bio[cols_bio].drop_duplicates(), how="left", on="Biopsy_Id")

    # select samples based on values if any
    if len(selection)>0:
        for col, val, mod in selection.items():
            if mod=="regex":
                r = re.compile(r'.*(%s).*' % val)
                mask = df_cln[col].apply(lambda x: bool(r.match(x)))
            elif mode=="exact":
                if type(val)!=list:
                    mask = df_cln[col]==val
                else:
                    mask = df_cln[col].isin(val)
            df_cln = df_cln[mask].copy()


    # select tumor type
    if args.tumor_type != "All":
        df_cln = df_cln.loc[df_cln[col_tt].isin(args.tumor_type.split("_"))]

    # choose column of drugs
    col_drug = "%s_Before_Biopsy" % level
    if col_drug not in df_cln.columns:
        raise ValueError("-ERROR: the level '%s' is not supported" % args.level)

    # remove patients with no data about drugs received
    mask_no_data = df_cln[col_drug].isnull()
    df_cln = df_cln.loc[~mask_no_data].copy()

    # if Classes, split by target. Drugs having multiple targets will be counted multiple times
    if level == "Classes":
        df_cln[col_drug] = df_cln[col_drug].apply(split_targeted_therapy_targets)

    # build dataframe of interest
    cols_keep = ["Subject_Id", col_tt, "Histological_Type", "Primary_Site", "Survival_Status", "Survival_Time",
                 col_drug] + cols_bio
    df_drug = df_cln[cols_keep].copy()

    # split aggregated drugs
    df_drug = explode_df(df=df_drug, cols=[col_drug], sep="|")

    return df_drug


def draw_double_heatmap_drugs(df_drug, col_drug):
    df_drug_sub = df_drug[["Subject_Id", col_drug]].copy()
    df_drug_sub["Received"] = 1
    df_drug_sub = df_drug_sub.loc[df_drug_sub[col_drug]!="None"]
    df_values = df_drug_sub.pivot(index="Subject_Id", columns=col_drug, values="Received").fillna(0)

    # organize the columns by hierarchical clustering using internal code from seaborn
    try:
        sns_res = sns.clustermap(df_values, metric="euclidean", standard_scale=1, method="ward")
        cols_ordered = sns_res.data2d.columns.tolist()
    except:
        cols_ordered = df_values.columns.tolist()
    df_values = df_values[cols_ordered]

    dfs = build_double_heatmap(df_values=df_values)

    brown_to_green_colors = ["#8c5322", "#c2842a", "#debf7c", "#f5ebc4", "#c8e7e5", "#7fcdc2", "#339a92",  "#11675f"]
    brown_to_green_cmap = cm.colors.LinearSegmentedColormap.from_list("BrownToGreen", colors=brown_to_green_colors)

    # get symbol sizes/fontsizes
    n_drugs = df_values.shape[1]
    fwer_size = max(2.5, 80 - 2*n_drugs + 0.0135*n_drugs**2)
    fdr_size = max(1, 11 - 0.15*n_drugs)
    ticks_labelsize = max(4, 11 - 0.1*n_drugs)

    config = DoubleHeatmapConfig()
    config.figure["figsize"] = (8, 8)
    config.heatmap["orientation"] = "antidiagonal"
    config.heatmap["orientation"] = "antidiagonal"
    config.heatmap["ticks_labelsize"] = ticks_labelsize
    config.count["cbar_xy"] = (0.01, 0.5)
    config.count["cbar_fraction"] = 0.35
    config.count["cbar_title"] = "Counts"
    config.ratio["cbar_fraction"] = 0.35
    config.ratio["cbar_title"] = "OR"
    config.ratio["cbar_xy"] = (0.5, 0.05)
    config.ratio["boundaries"] = [0, 0.001, 0.01, 0.1, 1, 5, 20, 100, 1000]
    config.ratio["cmap"] = brown_to_green_cmap
    config.test["fwer_size"] = fwer_size
    config.test["fdr_size"] = fdr_size

    fig, axes = plot_double_heatmap(df_count=dfs["count"],
                                    df_ratio=dfs["ratio"],
                                    df_test=dfs["test"],
                                    config=config)

    return fig


def main(args):
    col_drug = "%s_Before_Biopsy" % args.level
    col_tt = "Project_TCGA_More"

    # get table of drug classes
    df_drug_table = pd.read_excel(args.drug_table)

    # get data
    df_drug = get_table_of_drug_data(cohort=args.cohort, tumor_type=args.tumor_type, level=args.level)

    # special case: for tumor_type=All, consider only drugs given to at least 5 different patients
    if args.tumor_type=="All":
        df_cnt = df_drug[col_drug].value_counts()
        drugs_keep = df_cnt.loc[df_cnt>=5].index.tolist()
        drugs_drop = df_cnt.loc[df_cnt<5].index.tolist()
        df_drug = df_drug.loc[df_drug[col_drug].isin(drugs_keep)]
        if len(drugs_drop)>0:
            if len(drugs_drop)<20:
                print("-INFO: dropped following rare drugs:")
                print("\t"+ "\n\t".join(drugs_drop))
            else:
                print("-INFO: dropped %s rare drugs" % len(drugs_drop))

    # add classes if applicable
    if args.level == "Drugs":
        cols_lvl = [x for x in df_drug_table if "_Lvl_" in x]
        df_drug_table_unique = df_drug_table[["DCI"]+cols_lvl].drop_duplicates().rename(columns={"DCI": col_drug})
        df_drug = df_drug.merge(df_drug_table_unique, how="left", on=col_drug)

    # draw double heatmap showing cooccurrence of drugs
    fig = draw_double_heatmap_drugs(df_drug, col_drug)

    if args.level == "Classes":
        plt.subplots_adjust(bottom=0.35, left=0.35)
    else:
        plt.subplots_adjust(bottom=0.25, left=0.25)

    fig.savefig(args.output, dpi=300)
    print("-file saved at %s" % args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw barplots showing drug distribution.")
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--drug_table', type=str, help='Path to table of drugs.',
                        default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
    parser.add_argument('--tumor_type', type=str, default="PAAD", help='Selection of a subset of the cohort.')
    parser.add_argument('--level', type=str, default="Classes",
                        help='Chooose between "Classes" and "Drugs"')
    parser.add_argument('--output', type=str, help='Paths to output double heatmap plot.',
                        default="../../../results/treatment_resistances/plots/double_htmp_PAAD_Classes.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
