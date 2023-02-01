# -*- coding: utf-8 -*-
"""
@created: Feb 01 2022
@modified: Nov 06 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Intersect a table of CNV segments with a BED file of gene coordinates and add a filter status.
"""

import argparse
import gzip
import numpy as np
import os
import pandas as pd
import subprocess
import re

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
    df = pd.read_table(path, skiprows=len(header), na_values=["-","."])
    return df


def convert_num_to_str(x):
    try:
        y = "%d" % int(x)
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = x
        except:
            y = x

    return y


def main(args):
    if args.output.endswith(".tsv.gz"):
        pattern = ".tsv.gz"
    else:
        pattern = ".tsv"

    bed_b = args.output.replace(pattern, "_b.bed")
    bed_i = args.output.replace(pattern, "_i.bed")

    # load table and genes
    df_tab = read_table(args.input_tab)

    # make bed from tab
    cols_bed = ["chrom", "start", "end", "tcn.em", "lcn.em", "cf.em", "svtype", "svlen", "copy_number"] + \
            ["copy_number_more", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    df_bed = df_tab[cols_bed].copy()
    for col in ["tcn.em", "lcn.em"]:
        df_bed[col] = df_bed[col].apply(lambda x: "%d" % x if not np.isnan(x) else x)

    df_bed.fillna(".").to_csv(bed_b, index=False, header=False, sep="\t")

    # run bedtools intersect
    cmd = "bedtools intersect -a %s -b %s -wao > %s" % (args.input_bed, bed_b, bed_i)
    print("-running the command:\n\t%s" % cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # read results from bedtools intersect
    df_bed_i = pd.read_table(bed_i, header=None, sep="\t")
    df_bed_i.columns = ["chrom_gene", "start_gene", "end_gene", "gene_id", "gene_name", "gene_biotype", "gene_source"] \
            + cols_bed + ["overlap"]

    # set "." to NA and set everything to NA for genes that could not be intersected
    df_bed_i = df_bed_i.replace(".", np.nan)
    mask_genes_mis = df_bed_i["chrom"].isnull()
    df_bed_i.loc[mask_genes_mis, list(df_bed.columns)] = np.nan

    cols_old2new = {"start": "svstart", "end": "svend"}
    df_bed_i = df_bed_i.rename(columns=cols_old2new)

    # select columns keep
    cols_keep = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "chrom_gene", "start_gene", "end_gene",
                 "gene_id", "gene_name", "gene_biotype", "gene_source", "tcn.em", "lcn.em", "cf.em", "overlap",
                 "svtype", "svstart", "svend", "svlen", "copy_number", "copy_number_more"]
    cols_old2new = {"chrom_gene": "chrom", "start_gene": "start", "end_gene": "end", "gene_name": "gene"}

    df_cnv = df_bed_i[cols_keep]
    df_cnv = df_cnv.rename(columns=cols_old2new)

    # drop genes with overlap 0
    mask = df_cnv["overlap"]!=0
    df_cnv = df_cnv.loc[mask]
    print("-INFO: dropped %d/%d lines (~ genes) with 0 overlap" % (sum(~mask), len(mask)))

    # for genes with different copy-number, select the copy-number from the smallest SV (prioritize focal events)
    df_cnv = df_cnv.sort_values(by=["gene", "svlen"], ascending=True)
    n_row_bef = df_cnv.shape[0]
    df_cnv = df_cnv.drop_duplicates(subset=["gene"], keep="first")
    n_row_aft = df_cnv.shape[0]
    print("-INFO: dropped %d/%d lines from genes overlapping multiple SV" % (n_row_bef-n_row_aft, n_row_bef))

    df_cnv = df_cnv.rename(columns={"gene": "Hugo_Symbol", "gene_id": "Ensembl_Gene_Id",
                                    "chrom": "Chromosome", "copy_number": "Copy_Number",
                                    "copy_number_more": "Copy_Number_More"})
    df_cnv["TCN_EM:LCN_EM"] = df_cnv[["tcn.em", "lcn.em"]].astype(str).apply(":".join, axis=1)

    # add filter for events covering more than X Mb
    df_cnv[["svstart", "svend", "svlen"]]

    df_cnv["svlen"] = df_cnv["svlen"].apply(lambda x: re.sub(r"\.0+$", "", x))
    df_cnv["svlen"] = df_cnv["svlen"].astype(int)
    mask = df_cnv["svlen"] < args.threshold*1e6
    print("-INFO: flagged %d/%d lines (~ genes) from SV longer than %s Mb" % (sum(~mask), len(mask), args.threshold))
    df_cnv["FILTER"] = "PASS"
    if sum(~mask)>0:
        df_cnv.loc[~mask, "FILTER"] = "SV > %d Mb" % args.threshold

    # format numeric values
    cols_num = ["start", "end", "tcn.em", "lcn.em", "cf.em", "overlap", "svstart", "svend", "svlen", "Copy_Number"]
    for col_num in cols_num:
        if col_num in df_cnv:
            df_cnv[col_num] = df_cnv[col_num].fillna("NA").apply(convert_num_to_str)
    df_cnv["Copy_Number_More"] = df_cnv["Copy_Number_More"].fillna("NA")

    cols_gby = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Hugo_Symbol", "Ensembl_Gene_Id", "Chromosome"]
    cols_agg = ["Copy_Number", "Copy_Number_More", "TCN_EM:LCN_EM", "cf.em", "svtype", "svlen", "svstart", "svend",
                "overlap", "FILTER"]
    dt_agg = {x: ";".join for x in cols_agg}
    for col_gby in cols_gby:
        df_cnv[col_gby] = df_cnv[col_gby].fillna("NA")

    if len(df_cnv) > 0:
        df_cnv = df_cnv.groupby(cols_gby).agg(dt_agg).reset_index()
    else:
        df_cnv = df_cnv[cols_gby + cols_agg]

    # save and remove temporary files
    df_cnv.to_csv(args.output, index=False, sep="\t")
    os.remove(bed_b)
    os.remove(bed_i)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Used a bed file of genes to infer gene-level CNAs.")
    parser.add_argument('--input_tab', type=str, help='Path to tsv file.')
    parser.add_argument('--input_bed', type=str, help='Path to bed file.')
    parser.add_argument('--threshold', type=int, help='CNV events from segments larger than x Mb are flagged.')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
