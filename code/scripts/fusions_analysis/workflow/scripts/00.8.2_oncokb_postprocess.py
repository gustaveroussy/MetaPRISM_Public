# -*- coding: utf-8 -*-
"""
@created: Feb 28 2022
@modified: Feb 28 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Postoprocess fusion tables annotated by oncokb-annotator
"""

import argparse
import gzip
import numpy as np
import os
import pandas as pd

# functions ============================================================================================================

def read_header(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as file:
            header = [x for x in file.readlines() if x.startswith("##")]
    else:
        with open(path, "r") as file:
            header = [x for x in file.readlines() if x.startswith("##")]
    return header


def read_table(path):
    header = read_header(path)
    df = pd.read_table(path, skiprows=len(header), na_values=["."], sep="\t")
    return df


def convert_type_to_str(x):
    try:
        y = "%d" % int(x)
    except:
        try:
            y = "%d" % float(x)
            if y=="nan":
                y = np.nan
        except:
            y = x
    return y


def main(args):
    # load table and genes
    df_fus = read_table(args.table_fus)
    df_sam = read_table(args.table_sam)
    df_bio = read_table(args.table_bio)
    df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["RNA_N", "RNA_T"])]

    # select sample in the bio table
    mask = df_sam["Sample_Id"].isin(df_bio["Sample_Id"])
    df_sam = df_sam.loc[mask]

    old2new = {"TUMOR_SAMPLE_BARCODE": "Sample_Id", "FUSION": "Fusion_Id", "GENE1": "Gene_1", "GENE2": "Gene_2"}

    # split per sample
    dfs_fus_okb = []
    sam_okb = []
    for filepath in args.oncokb:
        sam = os.path.basename(filepath).split(".tsv")[0]
        sam_okb.append(sam)
        df_fus_okb = pd.read_table(filepath)
        df_fus_okb = df_fus_okb.rename(columns=old2new)
        dfs_fus_okb.append(df_fus_okb)

    sam_not_okb = list(set(df_sam["Sample_Id"]).difference(set(sam_okb)))
    for sam in sam_not_okb:
        filepath = os.path.join(os.path.dirname(args.oncokb[0]), "%s.tsv" % sam)
        print("-WARNING: %s does not exist" % filepath)

    # concatenate per sample table
    df_fus_okb = pd.concat(dfs_fus_okb)

    if "MUTATION_EFFECT" in df_fus_okb:
        df_fus_okb = df_fus_okb.loc[df_fus_okb["MUTATION_EFFECT"]!="Unknown"].copy()
    else:
        df_fus_okb = df_fus_okb.iloc[:0,:]

    if df_fus_okb.shape[0] > 0:
        df_fus_okb["Keep"] = "Y"
        cols_com = list(set(df_fus.columns).intersection(set(df_fus_okb.columns)))
        n_row_bef = df_fus.shape[0]
        df_fus = df_fus.merge(df_fus_okb, how="left", on=cols_com)
        n_row_aft = df_fus.shape[0]

        assert n_row_bef==n_row_aft

        # filter to fusions with non-empty annotations
        df_fus = df_fus.loc[df_fus["Keep"]=="Y"].copy()
        del df_fus["Keep"]

    # format some numeric columns
    cols = ["Breakpoint_1", "Breakpoint_2"] + [x for x in df_fus if x.endswith("_Reads")]
    for col in cols:
        df_fus[col] = df_fus[col].apply(convert_type_to_str)

    # save
    df_fus.to_csv(args.output, index=False, sep="\t")
    print("-file saved at %s" % args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare fusions for annotations by oncokb-annotator.")
    parser.add_argument('--table_fus', type=str, help='Path to fusion table.',
                        default="../../../data/prism/rna/fusions/prism_annotated_filtered.tsv.gz")
    parser.add_argument('--table_sam', type=str, help='Path to sample table.',
                        default="../../../data/prism/rna/fusions/sample_list.tsv")
    parser.add_argument('--table_bio', type=str, help='Path to biospecimen table.',
                        default="../../../data/prism/clinical/curated/bio_prism_in_design_curated.tsv")
    parser.add_argument('--oncokb', nargs="+", type=str, help='Paths to fusions annotated tables.')
    parser.add_argument('--output', type=str, help='Path to output table.',
                        default="../../../data/prism/rna/fusions/prism_annotated_filtered_oncokb.tsv.gz")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
