# -*- coding: utf-8 -*-
"""
@created: Jan 24 2022
@modified: Oct 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Clean oncokb annotations.
"""

import argparse
import numpy as np
import os
import pandas as pd

# functions ============================================================================================================

def main(args):
    # load table and rules
    df_table = pd.read_table(args.input)
    df_rules = pd.read_excel(args.rules, sheet_name="Oncokb_Rules")

    # build col merging MUTATION_EFFECT and ONCOGENIC
    cols = ["MUTATION_EFFECT", "ONCOGENIC"]
    df_table["Match"] = df_table[cols].fillna("N/A").apply("_".join, axis=1)
    df_rules["Match"] = df_rules[cols].fillna("N/A").apply("_".join, axis=1)

    # select
    df_table = df_table.merge(df_rules[["Match", "Keep"]], how="left", on="Match")
    mask = df_table["Keep"]=="Y"
    if (sum(~mask)>0):
        df_c = df_table["Match"].loc[~mask].value_counts()
        for name, count in df_c.items():
            print("-INFO: dropped %d rows with MUTATION_EFFECT_ONCOGENIC = %s" % (count, name))

    df_table = df_table.loc[df_table["Keep"]=="Y"].copy()
    del df_table["Match"]
    del df_table["Keep"]

    # split SAMPLE_ID if it contains aggregated tumor and normal sample ids
    if "SAMPLE_ID" in df_table:
        if df_table.shape[0]==0:
            df_table["Tumor_Sample_Barcode"] = np.nan
            df_table["Matched_Norm_Sample_Barcode"] = np.nan
            del df_table["SAMPLE_ID"]

        elif "_vs_" in df_table["SAMPLE_ID"].iloc[0]:
            cols = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
            df_table[cols] = df_table["SAMPLE_ID"].apply(lambda x: x.split("_vs_")).apply(pd.Series)
            del df_table["SAMPLE_ID"]

    if args.category == "mut":
        # OncoKB is likely misannotating variants from some categories. Consequently, only some
        # categories of variants are retained.
        if "Variant_Classification" not in df_table:
            print("-WARNING: Variant_Classification column is absent which prevents us from filtering potential" \
                  + " false-positives")
        else:
            vcs_keep = ["Missense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                        "Nonsense_Mutation", "Splice_Site", "Translation_Start_Site"]
            mask = df_table["Variant_Classification"].isin(vcs_keep)
            print("-INFO: selected %d/%d rows with Variant_Classification in:" % (sum(mask), len(mask)))
            print("\t" + "\n\t".join(vcs_keep))
            df_table = df_table.loc[mask].copy()

    # save
    df_table.to_csv(args.output, index=False, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Postprocess OncoKB anntoations.")
    parser.add_argument('--input', type=str, help='Path to input table.')
    parser.add_argument('--rules', type=str, help='Path to table of rules for cleaning.',
                        default="../../data/resources/oncokb/OncoKB_Curation_And_Rules.xlsx")
    parser.add_argument('--category', type=str, help='Choose one of cna, mut or fus.')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
