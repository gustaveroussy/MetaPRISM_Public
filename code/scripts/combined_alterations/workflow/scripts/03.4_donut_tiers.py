# -*- coding: utf-8 -*-
"""
@created: Jun 12 2022
@modified: Jun 12 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Draw donut plots describing fraction of samples with actionable/resistance alterations.
"""

import argparse
import os
import pandas as pd
import sys
sys.path.append("../../data_overview/workflow/functions")

from util import get_table_counts, draw_plot_donut_count

# functions ============================================================================================================

def main(args):
    # load data
    df_cln = pd.read_table(args.cln)
    df_alt = pd.read_table(args.alt)

    # select samples
    mask_dna = df_cln["Sample_Type"].apply(lambda x: "DNA_T" in x)
    mask_rna = df_cln["Sample_Type"].apply(lambda x: "RNA_T" in x)
    df_cln = df_cln.loc[mask_dna & mask_rna].copy()
    samples = df_cln["Sample_Id_DNA_T"].dropna().tolist() + df_cln["Sample_Id_RNA_T"].dropna().tolist()
    df_alt = df_alt.loc[df_alt["Sample_Id"].isin(samples)].copy()

    # add Tier status to each patient
    col_sub = "Subject_Id"
    col_tier = "%s_Level_Simple" % args.direction
    df_alt[col_tier] = df_alt[col_tier].fillna("Tier4")
    df_tiers = df_alt[[col_sub, col_tier]].sort_values(by=[col_sub, col_tier]).drop_duplicates(subset=col_sub)

    # merge with cln
    df_cln = df_cln.merge(df_tiers, how="left", on=col_sub)
    df_cln[col_tier] = df_cln[col_tier].fillna("Unknown").replace("Tier4", "Unknown")


    # table for plot
    df_count = df_cln[col_tier].value_counts().to_frame("Count").reset_index()
    df_count = df_count.rename(columns={"index": "Label"})

    # colors for plot
    if args.direction=="Sen":
        levels2colors = dict(Tier1="#00af56",
                             Tier2="#6fbfdd",
                             Tier3="#d4a8d0",
                             Unknown="#dadada")
    elif args.direction=="Res":
        levels2colors = dict(Tier1="#DD2D4A",
                             Tier2="#ffb700",
                             Tier3="#e0aaff",
                             Unknown="#dadada")

    df_count["Color"] = df_count["Label"].map(levels2colors)

    # draw and save
    fig = draw_plot_donut_count("PRISM", df_count, size=8, size_annotation=12)
    fig.write_image(args.output, width=125, height=125, engine="kaleido")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform a detail analysis of KRAS mutations in PAAD of TCGA.")
    parser.add_argument('--cln', type=str, help='Path to clinical file.',
                        default="../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
    parser.add_argument('--alt', type=str, help='Path to annotated alterations.',
                        default="../../../results/combined_alterations/alterations/aggregated_alterations_prism_all.tsv")
    parser.add_argument('--direction', type=str, help='Choose "Sen" or "Res".', default="Res")
    parser.add_argument('--output', type=str, help='Path to output file.',
                        default="../../../results/combined_alterations/donuts/all_samples_best_level.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
