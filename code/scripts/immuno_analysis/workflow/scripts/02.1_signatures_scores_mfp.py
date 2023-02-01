# -*- coding: utf-8 -*-
"""
@created: 26 Jul 2021
@modified: 15 Feb 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Learns a matrix of signatures scores from a matrix of gene expression.
"""

import argparse
import os
import numpy as np
import pandas as pd
from pathlib import Path
from pyprism.data import load_bio
import sys

filepath_mfp = "external/MFP/portraits"
sys.path.append(filepath_mfp)

from utils import read_gene_sets, ssgsea_formula

# functions ============================================================================================================

def remove_duplicates_summary(df_sum, bio_data):
    has_duplicates = df_sum["Sample_Id"].nunique() < df_sum.shape[0]
    is_tcga = all(df_sum["Col_Name"].str.startswith("TCGA"))
    if has_duplicates:
        if is_tcga:
            df_bio = pd.read_table(bio_data)
            df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["RNA_T", "RNA_N"])]
            cols_rna = [x for x in df_bio if x.startswith("RNA")]
            df_bio["RNA_All"] = df_bio[cols_rna].sum(axis=1)
            df_sum = df_sum.merge(df_bio[["Aliquot_Id", "RNA_All"]], how="left", on="Aliquot_Id")
            df_sum = df_sum.sort_values(by=["Sample_Id", "RNA_All"], ascending=False)
            df_sum = df_sum.drop_duplicates(subset=["Sample_Id"], keep="first")
            del df_sum["RNA_All"]
        else:
            print("-ERROR: unexpected duplicated ids in summary table")
    return df_sum


def load_expression(expression_data, expression_summary, bio_data):
    df_rna = pd.read_table(expression_data)
    df_sum = pd.read_table(expression_summary)
    df_sum = remove_duplicates_summary(df_sum, bio_data)
    cols_keep = ["ensembl_gene_id", "geneName", "geneID"] + df_sum["Col_Name"].tolist()
    df_rna = df_rna[[x for x in cols_keep if x in df_rna]]
    old2new = {r["Col_Name"]: r["Sample_Id"] for _, r in df_sum.iterrows()}
    df_rna = df_rna.rename(columns=old2new)

    return df_rna


def connect_gene_ids_to_symbols(df_rna, df_genes):
    df_rna = df_rna.rename(columns={"ensembl_gene_id": "Ensembl_Gene_Id_Version", "geneID": "Ensembl_Gene_Id_Version"})
    df_genes = df_genes[["Ensembl_Gene_Id_Version", "Hugo_Symbol"]]
    df_rna = df_rna.merge(df_genes, how="left", on="Ensembl_Gene_Id_Version")

    return df_rna


def select_genes(df_rna, df_genes, genesets, on="Ensembl_Gene_Id_Version"):
    set_symbols = list(set().union(*[geneset.genes for geneset in genesets.values()]))

    if on=="Ensembl_Gene_Id":
        df_rna["Ensembl_Gene_Id"] = df_rna["Ensembl_Gene_Id_Version"].apply(lambda x: x.split(".")[0])


    mask_id = df_rna[on].isin(df_genes[on])
    mask_symbol = df_rna["Hugo_Symbol"].isin(set_symbols)

    return df_rna.loc[mask_id | mask_symbol]


def clear_gene_ids(df_rna):
    for col in ["Ensembl_Gene_Id_Version", "Ensembl_Gene_Id", "geneID", "geneName"]:
        if col in df_rna.columns:
            del df_rna[col]

    return df_rna


def check_genesets_symbols(df_rna, genesets):
    set_symbols = list(set().union(*[geneset.genes for geneset in genesets.values()]))
    rna_symbols = df_rna["Hugo_Symbol"].unique().tolist()
    set_not_rna = set(set_symbols).difference(set(rna_symbols))

    if len(set_not_rna) > 0:
        print("-WARNING: gene symbols of some inputs signatures are not in the input expression table", flush=True)
        print("\t" + "\n\t".join(list(set_not_rna)), flush=True)
    else:
        print("-INFO: all gene symbols from gene sets found in expression table")


def summarize_duplicated_symbols(df_rna):
    if df_rna["Hugo_Symbol"].nunique() < df_rna.shape[0]:
        return df_rna.groupby("Hugo_Symbol").mean().T
    else:
        return df_rna.set_index("Hugo_Symbol").T


def main(args):
    df_rna = load_expression(args.expression_data, args.expression_summary, args.bio_data)
    genesets = read_gene_sets(args.signatures)
    df_genes_table = pd.read_table(args.gene_table)
    df_genes_selection = pd.read_table(args.gene_selection)

    df_rna = connect_gene_ids_to_symbols(df_rna, df_genes_table)
    df_rna = select_genes(df_rna, df_genes_selection, genesets=genesets, on="Ensembl_Gene_Id_Version")
    check_genesets_symbols(df_rna, genesets)
    df_rna = clear_gene_ids(df_rna)
    df_rna = summarize_duplicated_symbols(df_rna)

    df_scores = ssgsea_formula(df_rna, genesets)

    # save
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    df_scores.T.to_csv(args.output, index=True, sep="\t")
    print("-file saved at %s" % args.output)


# run ==============================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute signature scores using MFP ssgsea_formula function.')
    parser.add_argument('--expression_data', type=str, help="Path to input expression data table.",
                        default="../../../data/tcga/rna/kallisto-tximport/quantification_genes_1037_tpm_full.txt.gz")
    parser.add_argument('--expression_summary', type=str, help="Path to input expression summary table.",
                        default="../../../data/tcga/rna/summary/kallisto-tximport_quantification_genes_tpm_full.tsv")
    parser.add_argument('--bio_data', type=str, help="Path to biospecimen data table.",
                        default="../../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv")
    parser.add_argument('--gene_table', type=str, help="Path to gene table linking gene ids to gene symbols.",
                        default="../../../results/immuno_analysis/gene_tables/gencode_v27_updated.tsv")
    parser.add_argument('--gene_selection', type=str, help="Path to gene selection table.",
                default="../../../results/immuno_analysis/gene_tables/gencode_v27_coding_no-mt_no-histones_updated.tsv")
    parser.add_argument('--signatures', type=str, default="resources/bagaev_2021/signatures/gene_signatures.gmt",
                        help='Path to gene sets file in GMT format.')
    parser.add_argument('--output', type=str, help='Path to output table of signature scores.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    main(args)
