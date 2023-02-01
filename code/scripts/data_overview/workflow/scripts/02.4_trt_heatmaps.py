# -*- coding: utf-8 -*-
"""
@created: Oct 20 2021
@modified: Oct 20 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Heatmaps showing the distribution of treatments (drug names and classes) across tumor types.
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys

# local functions
sys.path.append("workflow/functions")
from load import load_cln_cohorts

# pyprism
from pyprism.data import load_colors, load_table
from pyprism.util import explode_df

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

# prettypy
from prettypy.heatmap import plot_heatmap, HeatmapConfig
from prettypy.heatmap import build_double_heatmap, plot_double_heatmap, DoubleHeatmapConfig

def main(args):
    # load cln table
    df_cln = load_cln_cohorts([args.cohort])
    df_cln = df_cln.rename(columns={"Project_TCGA_More": "Tumor_Type"})

    # filter out patients with no treatments data
    mask_rm = df_cln["Systematic_Treatment_Before_Biopsy"].isnull()
    # mask_na = df_cln["Drugs_Before_Biopsy"].isnull()
    # mask_com = df_cln["Comment_Drugs_Before_Biopsy"].isin(["Pas de dossier", "Pas de traitement prébiopsie"])
    # mask_rm = mask_na & mask_com
    df_cln = df_cln.loc[~mask_rm]

    # select tumor types
    df_count = df_cln["Tumor_Type"].value_counts()
    tumor_types = df_count.loc[df_count >= args.min_count_cols].index
    df_cln = df_cln.loc[df_cln["Tumor_Type"].isin(tumor_types)].copy()

    # load drugs table
    df_drug = load_table(args.drug_table)
    cols_levels = sorted([x for x in df_drug if "Class_Lvl" in x])
    df_drug[cols_levels] = df_drug[cols_levels].fillna("NA")
    df_drug["Class"] = df_drug[cols_levels].apply(" - ".join, axis=1)
    df_drug["Class"] = df_drug["Class"].apply(lambda x: x.replace(" - NA", ""))
    df_drug["Class"] = df_drug["Class"].replace(["", "NA"], np.nan)

    # delineate Drugs_Before_Biopsy 
    df_cln["Drugs_Before_Biopsy"] = df_cln["Drugs_Before_Biopsy"].fillna("No treatment")
    df_cln_drug = explode_df(df_cln, cols=["Drugs_Before_Biopsy"], sep="|")
    df_cln_drug["Drugs_Before_Biopsy"] = df_cln_drug["Drugs_Before_Biopsy"].replace("No treatment", np.nan)

    # add drug classes
    cols_drug = ["DCI", "Class"]
    df_cln_drug = df_cln_drug.merge(df_drug[cols_drug].drop_duplicates(),
                                    how="left", left_on="Drugs_Before_Biopsy", right_on="DCI")
    df_cln_drug = df_cln_drug.rename(columns={"Class": "Classes_Before_Biopsy"})

    ## plots
    cols = ["Drugs_Before_Biopsy", "Classes_Before_Biopsy"]
    titles = ["Drug and Tumor Type - Before Biopsy", "Drug class and Tumor_Type - Before Biopsy"]
    for col, title, width, height, output in zip(cols, titles, args.widths, args.heights, args.outputs):
        # prepare table for plot
        df_plot = df_cln_drug.groupby(["Tumor_Type"])[col].value_counts()
        df_plot = df_plot.unstack(level=-1).fillna(0).T
        df_plot = df_plot.loc[df_plot.sum(axis=1) > args.min_count_rows, :]
        df_plot = df_plot/df_plot.sum(axis=0)

        config = HeatmapConfig()
        config.figure["figsize"] = (width, height)
        config.figure["n_grid"] = 14
        config.heatmap["ticks_labelsize"] = 8
        config.heatmap["xticks_labelrotation"] = 90
        config.heatmap["yticklabels"] = True
        config.cbar["ticks_labelsize"] = 8
        config.cbar["title_fontsize"] = 10
        config.cbar["title"] = "Percent"
        config.cbar["title_pad"] = 10
        config.cbar["xy"] = (0,0.375)
        config.cbar["fraction"] = 0.25
        config.cbar["cmap"] = sns.color_palette("Reds", n_colors=3, as_cmap=True)
        config.cbar["boundaries"] = [0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1]

        fig, axes = plot_heatmap(df_plot, config)
        plt.title(title, fontsize=12, loc="left", fontweight="bold", pad=20)
        fig.savefig(output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Heatmaps showing treatments distributions across tumor types.')
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--min_count_cols', type=int, default=10,
                        help='Each column (i.e tumor type) with at least --min_count_cols tumors will be in the plot.')
    parser.add_argument('--min_count_rows', type=int, default=5,
                        help='Each row (i.e drug or class) with at least --min_count_rows will be in the plot.')
    parser.add_argument('--drug_table', type=str, help='Path to the drug table.',
                        default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
    parser.add_argument('--widths', type=int, nargs="+", help='Width of the outputs in inches.', default=[8, 12])
    parser.add_argument('--heights', type=int, nargs="+", help='Height of the outputs in pixels.', default=[16, 10])
    parser.add_argument('--outputs', type=str, nargs="+", help='Paths to output plots.',
                    default=["../../../results/data_overview/treatments/heatmaps/heatmap_drug_and_tumor_type_prism.pdf",
                             "../../../results/data_overview/treatments/heatmaps/heatmap_class_and_tumor_type_prism.pdf"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
