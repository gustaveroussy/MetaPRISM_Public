# -*- coding: utf-8 -*-
"""
@created: Mar 30 2022
@modified: May 28 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Venn plot showing the overlap between RNA and DNA.
"""

import argparse
import pandas as pd

# pyprism
from pyprism.data import load_cln
from pyprism.util import explode_df

# prettypy
import matplotlib.pyplot as plt
from prettypy.venn import plot_venn, VennConfig

def main(args):
    # load cln table
    df_cln = load_cln(args.cohort)

    # delineate sample type
    df_cln = explode_df(df_cln, cols=["Sample_Type"], sep="|")
    df_cln = df_cln.loc[df_cln["Sample_Type"]!="DNA_N"]
    df_cln["Sample_Type"] = df_cln["Sample_Type"].replace({"DNA_T": "WES", "RNA_T": "RNAseq"})

    config = VennConfig(figsize=(4,4), alpha=1,
                        arrow_r=0.4,
                        offset_label=[0.4,0],
                        arrow_connection_rad=0.2,
                        arrow_shrinkA=0,
                        arrow_shrinkB=15,
                        arrow_color="black",
                        arrow_position="count",
                        colors=["#F6A38E", "#D03711"])
    fig, ax = plot_venn(df=df_cln, col_set="Sample_Type", col_identifier="Subject_Id", config=config)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.savefig(args.output_paper, dpi=300)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--output', type=str, help='Path to output plot.',
                        default="../../../results/data_overview/misc/venn_molecular.pdf")
    parser.add_argument('--output_paper', type=str, help='Path to output plot.',
                        default="../../../results/figures_paper/FE1b.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
