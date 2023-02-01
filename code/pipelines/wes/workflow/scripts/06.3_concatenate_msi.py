# -*- coding: utf-8 -*-
"""
@created: Jan 28 2022
@modified: Apr 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Concatenate multiple msi mantis results tables.
"""

import argparse
import numpy as np
import os
import pandas as pd

# functions ============================================================================================================

def read_and_add_ids(filepath):
    basename = os.path.basename(filepath)
    tsample = basename.split("_vs_")[0]
    nsample = basename.split("_vs_")[1].split(".tsv")[0]
    df = pd.read_table(filepath)
    cols = ["Difference"]
    df = df[cols]

    if df.shape[0]==0:
        row_na = {col: np.nan for col in cols}
        df = df.append(row_na, ignore_index=True)
    else:
        df = df.iloc[[-1],:]
    df.insert(0, "Tumor_Sample_Id", tsample)
    df.insert(1, "Normal_Sample_Id", nsample)
    df = df.rename(columns={"Difference": "Stepwise_Difference"})
    return df

def main(args):
    # load tables
    dfs = [read_and_add_ids(filepath) for filepath in args.input]

    # folder = "somatic_msi_mantis"
    # files = [x for x in os.listdir(folder) if x.endswith(".tsv") and "kmer_counts" not in x]
    # dfs = []
    # for i, file in enumerate(files):
    #     filepath = os.path.join(folder, file)
    #     df = read_and_add_ids(filepath)
    #     dfs.append(df)
    #     if (i+1)%(len(files)//100)==0:
    #         print("-processed %d/%d files" % (i+1, len(files)), flush=True)

    # concatenate
    df = pd.concat(dfs, axis=0)

    # add Stable/Unstable
    threshold = 0.40
    df.loc[df["Stepwise_Difference"]>=threshold, "Status"] = "Unstable"
    df.loc[df["Stepwise_Difference"]<threshold, "Status"] = "Stable"

    # save
    df.to_csv(args.output, index=False, sep="\t")


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple MSI tables.")
    parser.add_argument('--input', nargs="+", type=str, help='Paths to msi tables.')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
