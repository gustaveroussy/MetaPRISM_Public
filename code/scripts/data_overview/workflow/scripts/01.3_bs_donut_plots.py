# -*- coding: utf-8 -*-
"""
@created: Sep 27 2021
@modified: Oct 14 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Donut plots showing the distribution of tumor types in the cohorts.
"""

import argparse
import os
import pandas as pd
import sys

# local functions
sys.path.append("workflow/functions")
from load import load_cln_cohorts, load_bio_cohorts
from util import get_table_counts, draw_plot_donut_count

# pyprism
from pyprism.data import load_colors

# plotly
import plotly.io as pio
pio.templates.default = "plotly_white"

# functions ============================================================================================================

def prepare_counts_biopsy_site(dfs_cln, BS2colors, x="Biopsy_Site"):
    """ Prepare count tables for fancy donut plot.

    Parameters
    ----------
    dfs_cln : dict
        Dict of dataframes
    BS2colors : dict
        Dict of colors
    x: str, optioal
        Name of the tumor type column

    Returns
    -------
    dfs_counts: dict
        Dict of dataframes with columns "Pull", Label" and "Color for the pie plot.
    """
    # vals_x_not_tcga = set().union(*[set(df[x]) for cohort, df in dfs_cln.items() if cohort!="tcga"])
    # vals_x_tcga_only = set(dfs_cln["tcga"][x]).difference(vals_x_not_tcga)
    dfs_counts = {}

    for cohort, df_cln in dfs_cln.items():
        df = get_table_counts(df_cln, x, add_pct=True)
        df = df.rename(columns={x: "Label"})
        # if cohort != "tcga":
        #     df["Pull"] = df["Label"].apply(lambda y: 0.28 if "Not_TCGA" in y else 0)
        # else:
        #     df["Pull"] = df["Label"].apply(lambda y: 0.25 if y in vals_x_tcga_only else 0)
        df["Color"] = df["Label"].map(BS2colors)
        dfs_counts[cohort] = df

    return dfs_counts


def main(args):
    # data
    dfs_cln = load_cln_cohorts(args.cohorts)
    dfs_bio = load_bio_cohorts(args.cohorts)

    # add Biopsy_Site
    for cohort in args.cohorts:
        df_cln = dfs_cln[cohort]
        df_bio = dfs_bio[cohort]

        # sometimes the dna and rna sample were extracted from different blocs resulting in different biopsy ids
        # however the biopsy site should be the same
        df_cln_d = df_cln.loc[~df_cln["Sample_Id_DNA_T"].isnull()].copy()
        df_cln_r = df_cln.loc[df_cln["Sample_Id_DNA_T"].isnull()].copy()

        cols = ["Sample_Id", "Biopsy_Site"]
        df_bio = df_bio[cols].drop_duplicates()
        df_cln_d = df_cln_d.merge(df_bio, how="left", left_on="Sample_Id_DNA_T", right_on="Sample_Id")
        df_cln_r = df_cln_r.merge(df_bio, how="left", left_on="Sample_Id_RNA_T", right_on="Sample_Id")
        df_cln = pd.concat((df_cln_d, df_cln_r), axis=0)

        bs_short = {"Spinal cord, cranial nerves, and other parts of central nervous system": \
                    "Spinal cord and other parts of CNS",
                    "Hematopoietic and reticuloendothelial systems": "Blood, bone marrow, and spleen",
                    "Connective, subcutaneous and other soft tissues": "Connnective and other soft tissues"}
        df_cln["Biopsy_Site"] = df_cln["Biopsy_Site"].replace(bs_short)
        dfs_cln[cohort] = df_cln

    # colors
    BS2colors = load_colors(sheet="Biopsy_Site")

    # count tumor types
    dfs_counts_bs = prepare_counts_biopsy_site(dfs_cln, BS2colors)

    # plot
    for cohort, output in zip(args.cohorts, args.outputs):
        fig = draw_plot_donut_count(cohort, dfs_counts_bs[cohort], size=8, size_annotation=11, rotation=355)
        fig.write_image(output, width=args.width, height=args.height, engine="kaleido")
        print("-file saved at %s" % output)

        if cohort=="prism":
            output = args.outputs[-1]
            fig.write_image(output, width=args.width, height=args.height, engine="kaleido")
            print("-file saved at %s" % output)

# Run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohorts', type=str, nargs="+", help='Names of the cohorts.',
                        default=["prism", "met500", "tcga"])
    parser.add_argument('--width', type=int, help='Width of the plot in pixels', default=300)
    parser.add_argument('--height', type=int, help='Height of the plot in pixels', default=300)
    parser.add_argument('--outputs', type=str, nargs="+", help='Paths to output plots.',
                        default=["../../../results/data_overview/biopsy_sites/donut_plot_tumor_types_prism.pdf",
                        "../../../results/data_overview/biopsy_sites/donut_plot_tumor_types_met500.pdf",
                        "../../../results/data_overview/biopsy_sites/donut_plot_tumor_types_tcga.pdf",
                        "../../../results/figures_paper/F1a_right.pdf"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
