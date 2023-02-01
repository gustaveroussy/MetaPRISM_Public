# -*- coding: utf-8 -*-
"""
@created: Jan 11 2022
@modified: Oct 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Preprocess CNA/Mut/Fusion table for annotation by oncokb-annotator
"""

import argparse
import gzip
import numpy as np
import os
import pandas as pd
import subprocess

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
    df = pd.read_table(path, skiprows=len(header), na_values=["-","."], sep="\t")
    return df


def select_genes(df_alt, alt_gene_name, df_gen, gen_gene_name, table_gen):
    # checks
    if not alt_gene_name in df_alt:
        raise ValueError("-%s is not a column of %s" % (alt_gene_name, alt))

    if not gen_gene_name in df_gen:
        raise ValueError("-%s is not a column of %s" % (gen_gene_name, genes))

    # select genes
    if df_alt.shape[0] > 0:
        mask = df_alt[alt_gene_name].isin(df_gen[args.gen_gene_name])
        print("-INFO: selected %d/%d rows from genes in %s" % (sum(mask), len(mask), table_gen))

        return df_alt.loc[mask].copy()
    else:
        return df_alt


def load_and_preprocess_cna(table_alt, table_gen, table_cln, gen_gene_name):
    # load alterations
    df_cna = pd.concat([read_table(tab) for tab in table_alt])
    df_gen = read_table(table_gen)
    df_cln = read_table(table_cln)

    # add sampe pair
    col_pair = "DNA_P"
    col_tsb = "Tumor_Sample_Barcode"
    col_nsb = "Matched_Norm_Sample_Barcode"
    df_cna[col_pair] = df_cna[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
    print("-INFO: loaded %d lines from %d sample pairs" % (len(df_cna), len(table_alt)))

    # retain only calls with Copy_Number_More = -2 or 2
    mask_in = df_cna["Copy_Number"].isin([-2, 2])
    df_cna = df_cna.loc[mask_in]
    print("-INFO: selected %d/%d lines having Copy_Number at 2 or -2" % (sum(mask_in), len(mask_in)))
    cna_gene_name = "Hugo_Symbol"
    df_cna["Copy_Number"] = df_cna["Copy_Number"].astype(int)

    # select genes
    df_gen = df_gen.loc[df_gen["OncoKB Annotated"]=="Yes"]
    df_cna_okb = select_genes(df_cna, "Hugo_Symbol", df_gen, gen_gene_name, table_gen)

    # some symbols may correspond to two gene ids or more. in this case, select one at random
    cols = ["DNA_P", "Hugo_Symbol", "Copy_Number"]
    df_cna_okb = df_cna_okb.drop_duplicates(subset=cols, keep="first")

    # format for Oncokb
    df_cna_okb["Locus ID"] = 0
    df_cna_okb["Cytoband"] = 0
    cols_idx = ["Hugo_Symbol", "Locus ID", "Cytoband"]
    df_cna_okb = df_cna_okb.pivot(index=cols_idx, columns="DNA_P", values="Copy_Number")
    df_cna_okb = df_cna_okb.fillna(0).astype(int).reset_index(drop=False)

    df_cln_okb = df_cln[["DNA_P", "MSKCC_Oncotree"]]
    df_cln_okb = df_cln_okb.loc[df_cln_okb["DNA_P"].isin(df_cna_okb.columns[3:].tolist())]
    df_cln_okb = df_cln_okb.rename(columns={"DNA_P": "SAMPLE_ID", "MSKCC_Oncotree": "ONCOTREE_CODE"})

    return df_cna_okb, df_cln_okb


def load_and_preprocess_mut(table_alt, table_gen, table_cln, gen_gene_name, maf_ncbi_build):
    # load alterations
    df_mut = pd.concat([read_table(tab) for tab in table_alt])
    df_gen = read_table(table_gen)
    df_cln = read_table(table_cln)

    # info
    print("-INFO: loaded %d lines from %d sample pairs" % (len(df_mut), len(table_alt)))

    # add missing columns (this happens when the maf file is empty)
    cols_req = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele1",
                "Tumor_Seq_Allele2", "Hugo_Symbol"]

    for col in cols_req:
        if col not in df_mut:
            df_mut[col] = np.nan

    # cols fillna by "-"
    cols_fillna = ["Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"]
    df_mut[cols_fillna] = df_mut[cols_fillna].fillna("-")

    # set to str type
    convert_to_str = lambda x: "%d" % x if not np.isnan(x) else x
    df_mut["Start_Position"] = df_mut["Start_Position"].apply(convert_to_str)
    df_mut["End_Position"] = df_mut["End_Position"].apply(convert_to_str)

    # udpate contents of NCBI_Build
    df_mut["NCBI_Build"] = maf_ncbi_build

    # select genes
    df_gen = df_gen.loc[df_gen["OncoKB Annotated"]=="Yes"]
    df_mut_okb = select_genes(df_mut, "Hugo_Symbol", df_gen, gen_gene_name, table_gen)

    # format clinical table
    df_cln_okb = df_cln[["DNA_T", "MSKCC_Oncotree"]]
    df_cln_okb = df_cln_okb.rename(columns={"DNA_T": "SAMPLE_ID", "MSKCC_Oncotree": "ONCOTREE_CODE"})

    if df_mut_okb.shape[0] > 0:
        df_cln_okb = df_cln_okb.loc[df_cln_okb["SAMPLE_ID"].isin(df_mut_okb["Tumor_Sample_Barcode"].unique().tolist())]

        # add ONCOTREE_CODE to alterations table 
        df_cln_okb_m = df_cln_okb[["SAMPLE_ID", "ONCOTREE_CODE"]].drop_duplicates()
        df_cln_okb_m = df_cln_okb_m.rename(columns={"SAMPLE_ID": "Tumor_Sample_Barcode"})
        df_mut_okb = df_mut_okb.merge(df_cln_okb_m, how="left", on="Tumor_Sample_Barcode")
    else:
        df_cln_okb = df_cln_okb.iloc[:0,:]

    return df_mut_okb, df_cln_okb


def load_and_preprocess_fus(table_alt, table_gen, table_cln, gen_gene_name):
    # load alterations
    df_fus = pd.concat([read_table(tab) for tab in table_alt])
    df_gen = read_table(table_gen)
    df_cln = read_table(table_cln)

    # rename some columns
    df_fus_okb = df_fus.rename(columns={"Fusion_Id": "Fusion", "Sample_Id": "Tumor_Sample_Barcode", "Gene_1": "Gene1",
                                    "Gene_2": "Gene2"})

    # add missing columns (this happens when the maf file is empty)
    cols_req = ["Tumor_Sample_Barcode", "Fusion", "Gene1", "Gene2"]

    for col in cols_req:
        if col not in df_fus:
            df_fus_okb[col] = np.nan

    # format clinical table
    df_cln_okb = df_cln[["DNA_T", "MSKCC_Oncotree"]]
    df_cln_okb = df_cln_okb.loc[df_cln_okb["DNA_T"].isin(df_mut_okb["Tumor_Sample_Barcode"].unique().tolist())]
    df_cln_okb = df_cln_okb.rename(columns={"DNA_T": "SAMPLE_ID", "MSKCC_Oncotree": "ONCOTREE_CODE"})

    return df_fus_okb, df_cln_okb


def main(args):
    # load
    if args.category == "cna":
        df_alt, df_cln = load_and_preprocess_cna(table_alt=args.table_alt, table_gen=args.table_gen,
                                                 table_cln=args.table_cln, gen_gene_name=args.gen_gene_name)
    elif args.category == "mut":
        df_alt, df_cln = load_and_preprocess_mut(table_alt=args.table_alt, table_gen=args.table_gen,
                                                 table_cln=args.table_cln, gen_gene_name=args.gen_gene_name,
                                                 maf_ncbi_build=args.maf_ncbi_build)
    elif args.category == "fus":
        df_alt, df_cln = load_and_preprocess_fus(table_alt=args.table_alt, table_gen=args.table_gen,
                                                 table_cln=args.table_cln, gen_gene_name=args.gen_gene_name)
    else:
        raise ValueError("Unsupported value '%s' for --category. Choose 'cna', 'mut' or 'fus'" % args.category)

    # save
    df_alt.to_csv(args.output_alt, index=False, sep="\t")
    print("-table saved at %s" % args.output_alt)

    if args.output_cln is not None:
        df_cln.to_csv(args.output_cln, index=False, sep="\t")
        print("-table saved at %s" % args.output_cln)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset and concatenate CNV, MAF or fusion tables.")
    parser.add_argument('--table_alt', type=str, nargs="+", help='Path to alterations table(s).')
    parser.add_argument("--table_cln", type=str, help="Path to table of clinical data per each pair.")
    parser.add_argument('--table_gen', type=str, help='Path to table of genes.')
    parser.add_argument('--gen_gene_name', type=str, help='Column containing gene names in genes.',
                       default="Hugo Symbol")
    parser.add_argument('--category', type=str, help='Choose one of cna, mut or fus.')
    parser.add_argument('--maf_ncbi_build', type=str, default=None,
                        help='NCBI build for updating NCBI_Build value of MAF table.')
    parser.add_argument('--output_alt', type=str, help='Path to output table alterations table for oncokb annotator.')
    parser.add_argument('--output_cln', type=str, default=None, help='Path to output table clinical table for oncokb annotator.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
