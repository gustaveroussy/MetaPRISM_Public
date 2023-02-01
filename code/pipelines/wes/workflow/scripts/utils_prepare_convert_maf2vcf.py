# -*- coding: utf-8 -*-
"""
@created: Jan 05 2022
@modified: Oct 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Prepare MC3 MAF file or subsets of it to be converted into VCF via vcf2maf.
"""

import argparse
import numpy as np
import pandas as pd

def sort_maf(df):
    chr_all = df["Chromosome"].astype(str).unique().tolist()
    chr_normal = ["%d" % i for i in range(1,23)] + ["X", "Y"]
    chr_other = sorted([x for x in chr_all if x not in chr_normal])
    chr_old2new = {**{"X": "23", "Y": "24"}, **{c: "%d" % (i+25) for i, c in enumerate(chr_other)}}
    chr_new2old = {int(v):k for k,v in chr_old2new.items()}

    df["Chromosome"] = df["Chromosome"].astype(str).replace(chr_old2new).astype(int)
    df["Start_Position"] = df["Start_Position"].astype(float)
    df["End_Position"] = df["End_Position"].astype(float)
    df = df.sort_values(by=["Chromosome", "Start_Position", "Variant_Type"])
    df["Chromosome"] = df["Chromosome"].replace(chr_new2old).astype(str)
    df["Start_Position"] = df["Start_Position"].apply(lambda x: "%d" % x if not np.isnan(x) else x)
    df["End_Position"] = df["End_Position"].apply(lambda x: "%d" % x if not np.isnan(x) else x)
    return df


def change_sep_filter(df):
    df["FILTER"] = df["FILTER"].apply(lambda x: x.replace(",", ";"))
    return df


def caller_filtering_ind(x):
    if "INDELOCATOR" in x["CENTERS"] or "VARSCANI" in x["CENTERS"]:
        return x["FILTER"]
    else:
        if x["FILTER"]=="PASS":
            return "INDEL_CALLING_ALIGNMENT"
        else:
            return x["FILTER"] + ";" + "INDEL_CALLING_ALIGNMENT"


def caller_filtering_snp(x):
    if x["NCALLERS"] < 2:
        if x["FILTER"]=="PASS":
            return "SNP_CALLING_ALIGNMENT"
        else:
            return x["FILTER"] + ";" + "SNP_CALLING_ALIGNMENT"
    else:
        return x["FILTER"]


def add_caller_filtering(df):
    df["Order"] = np.arange(0, df.shape[0])
    df_ind = df.loc[df["Variant_Type"].isin(["INS", "DEL"])].copy()
    df_snp = df.loc[~df["Variant_Type"].isin(["INS", "DEL"])].copy()

    # for indels, filter out calls not found by INDELOCATOR or VARSCANI
    # sometimes, the caller name is suffixed by a *, that's not a problem
    df_ind["FILTER"] = df_ind[["FILTER", "CENTERS"]].apply(caller_filtering_ind, axis=1)

    # for snps, filter out calls found only by 1 caller.
    df_snp["FILTER"] = df_snp[["FILTER", "NCALLERS"]].apply(caller_filtering_snp, axis=1)

    df_new = pd.concat((df_ind, df_snp), axis=0)
    df_new = df_new.sort_values(by="Order")
    del df_new["Order"]

    return df_new


def drop_missing_alleles(df):
    mask_bad = (df["Tumor_Seq_Allele1"].isnull())|(df["Tumor_Seq_Allele1"]==".")
    mask_bad = mask_bad|(df["Tumor_Seq_Allele2"].isnull())|(df["Tumor_Seq_Allele2"]==".")
    mask_bad = mask_bad|(df["Reference_Allele"].isnull())|(df["Reference_Allele"]==".")
    if sum(mask_bad)>0:
        print("-warning! removing %d mutations with empty allele in ref or tumor" % sum(mask_bad))
        df = df.loc[~mask_bad].copy()
    return df


def main(args):
    df = pd.read_table(args.input, sep="\t")
    df = sort_maf(df)
    df = change_sep_filter(df)
    df = add_caller_filtering(df)
    df = drop_missing_alleles(df)
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correct line end delimited")
    parser.add_argument("--input", type=str, help="Path to input table")
    parser.add_argument("--output", type=str, help="Path to output table")

    args = parser.parse_args()
    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
