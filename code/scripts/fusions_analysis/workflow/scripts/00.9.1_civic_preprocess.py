# -*- coding: utf-8 -*-
"""
@created: Mar 01 2022
@modified: Mar 01 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Preprocess fusion table for annotation by CivicAnnotator
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


def main(args):
    # load table and genes
    df_fus = read_table(args.table_fus)
    df_sam = read_table(args.table_sam)
    df_bio = read_table(args.table_bio)
    df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["RNA_N", "RNA_T"])]

    # select sample in the bio table
    mask = df_sam["Sample_Id"].isin(df_bio["Sample_Id"])
    df_sam = df_sam.loc[mask]
    print("-INFO: %d/%d samples are selected as they are present in %s" % (sum(mask), mask.shape[0], args.table_bio))

    # rename some columns
    # add missing columns (this happens when the input file is empty)
    cols_req = ["Sample_Id", "Fusion_Id", "Gene_1", "Gene_2"]

    for col in cols_req:
        if col not in df_fus:
            df_fus[col] = np.nan

    # split per sample
    folder = args.output
    os.makedirs(folder, exist_ok=True)
    for sam in df_sam["Sample_Id"]:
        filepath = os.path.join(folder, "%s.tsv" % sam)
        df_fus_sam = df_fus.loc[df_fus["Sample_Id"]==sam, cols_req].drop_duplicates()
        df_fus_sam.to_csv(filepath, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare fusions for annotations by CivicAnnotator.")
    parser.add_argument('--table_fus', type=str, help='Path to fusion table.',
                        default="../../../data/prism/rna/fusions/prism_annotated_filtered.tsv.gz")
    parser.add_argument('--table_sam', type=str, help='Path to sample table.',
                        default="../../../data/prism/rna/fusions/sample_list.tsv")
    parser.add_argument('--table_bio', type=str, help='Path to biospecimen table.',
                        default="../../../data/prism/clinical/curated/bio_prism_in_design_curated.tsv")
    parser.add_argument('--output', type=str, help='Path to output folder.',
                        default="../../../data/prism/rna/fusions/civic_pre")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
