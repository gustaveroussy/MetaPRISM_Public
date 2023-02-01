# -*- coding: utf-8 -*-
"""
@created: 28 Jul 2021
@modified: 16 Feb 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Performs a selection of columns on the TCGA signature tables. Select only the samples for which a label immune subtypes
is available in the files provided by Bagaev et al.
"""

import argparse
import numpy as np
import os
import pandas as pd
from pathlib import Path

from pyprism.data import load_table, load_bio

# functions ============================================================================================================

def main(args):
    df_X = pd.read_table(args.input_signatures, index_col=0)
    df_Y = pd.read_table(args.input_annotation, index_col=0)
    df_bio = pd.read_table(args.input_bio)
    df_bio = df_bio.loc[df_bio["Sample_Type"]=="RNA_T"]
    df_bio = df_bio.loc[df_bio["Subject_Id"].isin(df_Y.index)]
    df_X = df_X[df_bio["Sample_Id"].unique()]

    # methods of the paper explicitly mention that patients with multiple rna-seq files were excluded
    # therefore, cases of patients with multiple samples are removed
    df_s = pd.DataFrame({"Sample_Id": df_X.columns})
    df_s["Subject_Id"] = df_s["Sample_Id"].str[:12]
    samples_to_be_removed  = df_s.loc[df_s["Subject_Id"].duplicated(keep=False)]["Sample_Id"].unique().tolist()
    df_X = df_X[[x for x in df_X if x not in samples_to_be_removed]]

    old_columns = df_X.columns
    new_columns = list(map(lambda x: x[:12], old_columns))

    df_cols = pd.DataFrame({"Old_Col": old_columns, "New_Col": new_columns})
    df_X = df_X[df_cols["Old_Col"]]
    df_Y = df_Y.loc[df_cols["New_Col"]]
    df_Y.index = df_cols["Old_Col"]
    df_Y.index.name = None

    # save
    os.makedirs(str(Path(args.output_signatures).parent.absolute()), exist_ok=True)
    df_X.to_csv(args.output_signatures, index=True, sep="\t")
    print("-file saved at %s" % args.output_signatures)
    df_Y.to_csv(args.output_annotation, index=True, sep="\t")
    print("-file saved at %s" % args.output_annotation)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Performs a selection of columns on the TCGA signature tables.')
    parser.add_argument('--input_signatures', type=str, help='Path to able of signatures scores per sample.',
                        default="../../../results/immuno_analysis/tcga/tables/signatures_bagaev_2021_mfp.tsv")
    parser.add_argument('--input_annotation', type=str,
                        help="Path to table of samples annotations containing immune subtypes in 'MFP' column.",
                        default="resources/bagaev_2021/Pan_TCGA/annotation.tsv")
    parser.add_argument('--input_bio', type=str, help="Path to table of samples data.",
                        default="../../../data/tcga/clinical/curated_other/bio_tcga_all_curated.tsv")
    parser.add_argument('--output_signatures', type=str, help='Output table of signatures scores per sample.',
                default="../../../results/immuno_analysis/tcga/tables/signatures_mfp_model.tsv")
    parser.add_argument('--output_annotation', type=str,
                default="../../../results/immuno_analysis/tcga/tables/annotation_mfp_model.tsv",
                help="Output table of samples annotations containing immune subtypes in 'MFP' column.")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
