# -*- coding: utf-8 -*-
"""
@created: Jan 28 2022
@modified: Nov 04 2022
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
import gzip
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


def read_and_add_ids(filepath):
    basename = os.path.basename(filepath)
    tsample = basename.split("_vs_")[0]
    nsample = basename.split("_vs_")[1].split(".vcf")[0]
    header = read_header(filepath)

    cols_old = ["purity", "ploidy", "dipLogR", "est_insert_size", "emflags"]
    cols_new = ["Purity", "Ploidy", "DipLogR", "Est_Insert_Size", "EM_Flags"]
    lines = {}

    for line in header:
        for col in cols_old:
            pattern = "##%s=" % col
            if line.startswith(pattern):
                try:
                    if col!="emflags":
                        value = float(line.split(pattern)[1])
                    else:
                        value = line.split(pattern)[1].strip()
                except:
                    value = np.nan

                lines[col] = [value]

    df = pd.DataFrame(lines)
    df = df.rename(columns={o:n for o,n in zip(cols_old, cols_new)})

    df.insert(0, "Tumor_Sample_Id", tsample)
    df.insert(1, "Normal_Sample_Id", nsample)
    return df


def main(args):
    # load tables
    dfs = [read_and_add_ids(filepath) for filepath in args.input]
    # folder = "results_prism/calling/somatic_cnv_facets"
    # files = [x for x in os.listdir(folder) if x.endswith(".vcf.gz")]
    # dfs = []
    # for i, file in enumerate(files):
    #     filepath = os.path.join(folder, file)
    #     df = read_and_add_ids(filepath)
    #     dfs.append(df)
    #     if (i+1)%(len(files)//100)==0:
    #         print("-processed %d/%d files" % (i+1, len(files)), flush=True)

    # concatenate
    df = pd.concat(dfs, axis=0)

    # save
    df.to_csv(args.output, index=False, sep="\t")


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate sample-level metrics computed by CNV facets.")
    parser.add_argument('--input', nargs="+", type=str, help='Paths to cnv facets vcf files.')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
