# -*- coding: utf-8 -*-
"""
@created: Jun 10 2022
@modified: Jun 10 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Donut plots showing capture kits used in WES and RNAseq.
"""

import argparse
import pandas as pd
import sys

# pyprism
from pyprism.data import load_cln, load_bio

# local functions
sys.path.append("workflow/functions")
from util import get_table_counts, draw_plot_donut_count


def main(args):
    # load cln table
    df_cln = load_cln(args.cohort)
    df_bio = load_bio(args.cohort)

    # split dna and rna sample
    df_cln_dna = df_cln.loc[df_cln["Sample_Type"].apply(lambda x: "DNA_T" in x)].copy()
    df_cln_dna["Sample_Id"] = df_cln_dna["Sample_Id_DNA_T"]
    df_cln_dna = df_cln_dna.drop_duplicates(subset=["Sample_Id"])
    df_cln_rna = df_cln.loc[df_cln["Sample_Type"].apply(lambda x: "RNA_T" in x)].copy()
    df_cln_rna["Sample_Id"] = df_cln_rna["Sample_Id_RNA_T"]
    df_cln_rna = df_cln_rna.drop_duplicates(subset=["Sample_Id"])

    # add bio Capture_Kit, Platform, Library_Selection
    cols_bio = ["Sample_Id", "Capture_Kit", "Platform", "Library_Selection"]
    df_cln_dna = df_cln_dna.merge(df_bio[cols_bio].drop_duplicates(), how="left", on="Sample_Id")
    df_cln_rna = df_cln_rna.merge(df_bio[cols_bio].drop_duplicates(), how="left", on="Sample_Id")


    # dna donut
    col_dna = "Capture_Kit"
    values2colors = dict(SureSelect_CR="#e4c1f9",
                         SureSelect_CR2="#ff99c8",
                         SureSelect_V5="#a9def9",
                         TWISBioscienceCustomIG_Targeted="#fcf6bd",
                         TWISBioscienceCustomIG_TargetedV2="#d0f4de")
    df_count = get_table_counts(df_cln_dna, col_dna, add_pct=True)
    df_count = df_count.rename(columns={col_dna: "Label"})
    df_count["Color"] = df_count["Label"].map(values2colors)

    fig = draw_plot_donut_count("WES", df_count, size=7, size_annotation=11, textposition="outside")
    fig.write_image(args.outputs[0], width=300, height=300, engine="kaleido")

    # rna donut
    col_rna = "Library_Selection"
    values2colors = dict(PolyA="#d0f4de")
    df_count = get_table_counts(df_cln_rna, col_rna, add_pct=True)
    df_count = df_count.rename(columns={col_rna: "Label"})
    df_count["Color"] = df_count["Label"].map(values2colors)

    fig = draw_plot_donut_count("RNAseq", df_count, size=15, size_annotation=24, textposition="outside")
    fig.write_image(args.outputs[1], width=300, height=300, engine="kaleido")


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--outputs', type=str, nargs="+", help='Path to output plots.',
                        default=["../../../results/data_overview/misc/donut_WES.pdf",
                                 "../../../results/data_overview/misc/donut_RNAseq.pdf"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
