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
from load import load_cln_cohorts
from util import get_table_counts, draw_plot_donut_count

# pyprism
from pyprism.data import load_colors

# plotly
import plotly.io as pio
import plotly.graph_objects as go
pio.templates.default = "plotly_white"

# functions ============================================================================================================

def prepare_counts_tumor_type(dfs_cln, PTM2colors, x="Project_TCGA_More"):
    """ Prepare count tables for fancy donut plot.

    Parameters
    ----------
    dfs_cln: dict
        Dict of dataframes
    PTM2colors:
        Dict of colors
    x: str, optioal
        Name of the tumor type column

    Returns
    -------
    dfs_counts: dict
        Dict of dataframes with columns "Pull", Label" and "Color for the pie plot.
    """
    vals_x_not_tcga = set().union(*[set(df[x]) for cohort, df in dfs_cln.items() if cohort!="tcga"])
    vals_x_tcga_only = set(dfs_cln["tcga"][x]).difference(vals_x_not_tcga)
    dfs_counts = {}

    for cohort, df_cln in dfs_cln.items():
        df = get_table_counts(df_cln, x, add_pct=True)
        df = df.rename(columns={x: "Label"})
        if cohort != "tcga":
            df["Pull"] = df["Label"].apply(lambda y: 0.28 if "Not_TCGA" in y else 0)
        else:
            df["Pull"] = df["Label"].apply(lambda y: 0.25 if y in vals_x_tcga_only or "Not_TCGA" in y else 0)
        df["Color"] = df["Label"].map(PTM2colors)
        dfs_counts[cohort] = df

    return dfs_counts


def main(args):
    # data
    dfs_cln = load_cln_cohorts(args.cohorts)
    for cohort, df_cln in dfs_cln.items():
        df_cln["Project_TCGA_More"] = df_cln["Project_TCGA_More"].replace({"Unknown_Primary": "Unknown primary"})

    # colors
    PTM2colors = load_colors(sheet="Project_TCGA_More")

    # count tumor types
    dfs_counts_tt = prepare_counts_tumor_type(dfs_cln, PTM2colors)

    # plot
    for cohort, output, output_paper in zip(args.cohorts, args.outputs, args.outputs_paper):
        if cohort=="tcga":
            width = 0.8*args.width
            height = 0.8*args.height
        else:
            width = args.width
            height = args.height
        fig = draw_plot_donut_count(cohort, dfs_counts_tt[cohort], size=8, size_annotation=11)
        fig.write_image(output, width=width, height=height, engine="kaleido")
        fig.write_image(output_paper, width=width, height=height, engine="kaleido")
        print("-file saved at %s" % output)
        print("-file saved at %s" % output_paper)

# Run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohorts', type=str, nargs="+", help='Names of the cohorts.',
                        default=["prism", "met500", "tcga"])
    parser.add_argument('--width', type=int, help='Width of the plot in pixels', default=300)
    parser.add_argument('--height', type=int, help='Height of the plot in pixels', default=300)
    parser.add_argument('--outputs', type=str, nargs="+", help='Paths to output plots.',
                        default=["../../../results/data_overview/donut_plot_tumor_types_prism.pdf",
                        "../../../results/data_overview/donut_plot_tumor_types_met500.pdf",
                        "../../../results/data_overview/donut_plot_tumor_types_tcga.pdf"])
    parser.add_argument('--outputs_paper', type=str, nargs="+", help='Paths to output plots.',
                        default=["../../../results/figures_paper/F1a.pdf",
                                 "../../../results/figures_paper/FS3c.pdf",
                                 "../../../results/figures_paper/FS3a.pdf"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
