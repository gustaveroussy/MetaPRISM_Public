# -*- coding: utf-8 -*-
"""
@created: Jul 27 2021
@modified: Feb 14 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Update gene symbols in table of genes.
"""

import argparse
import pandas as pd
from pyprism.data import load_table
from pyprism.util import explode_df

# function =============================================================================================================

def update_gene_symbol_using_hgnc(df_gene, df_hgnc, gene_name=None, gene_ensembl_id=None, chr_name=None, verbose=True):
    cols_gene = [x for x in df_gene]
    cols_hgnc = ["Approved symbol", "Ensembl gene ID"]
    df_hgnc_approved = df_hgnc.loc[df_hgnc["Status"]=="Approved"].copy()
    df_hgnc_approved["Chr"] = df_hgnc_approved["Chromosome"].str.extract("^([0-9A-Z]+)")

    cols_on = ["Ensembl gene ID", "Approved symbol", "Previous symbols", "Alias symbols"]
    dfs_gene = {}
    genes_updated = []

    for col_on in cols_on:
        if col_on in cols_hgnc:
            cols_merge = cols_hgnc
        else:
            cols_merge = cols_hgnc + [col_on]

        if col_on=="Ensembl gene ID":
            right_on = col_on
            left_on = gene_ensembl_id
        else:
            right_on = col_on
            left_on = gene_name

        df_hgnc_sub = df_hgnc_approved.loc[~df_hgnc_approved[col_on].isnull()]
        df_hgnc_sub = explode_df(df_hgnc_sub, cols=[col_on], sep=",")
        df_gene_sub = df_gene.loc[~df_gene[gene_name].isin(genes_updated)]

        if chr_name in df_gene:
            df_right = df_hgnc_sub[cols_merge+["Chr"]]
            df_right = df_right.drop_duplicates(subset=[col_on, "Chr"], keep="first")

            df_gene_sub_chr = df_gene_sub.loc[~df_gene_sub[chr_name].isnull()]
            df_gene_sub_nochr = df_gene_sub.loc[df_gene_sub[chr_name].isnull()]

            df_merge_chr = df_gene_sub_chr.merge(df_right, how="inner", left_on=[left_on, chr_name],
                                                 right_on=[col_on, "Chr"])
            df_merge_no_chr = df_gene_sub_nochr.merge(df_right, how="inner", left_on=left_on, right_on=col_on)
            df_merge = pd.concat((df_merge_chr, df_merge_no_chr), axis=0)
        else:
            df_right = df_hgnc_sub.drop_duplicates(subset=[col_on], keep="first")
            df_merge = df_gene_sub.merge(df_right, how="inner", left_on=left_on, right_on=col_on)

        dfs_gene[col_on] = df_merge
        genes_updated = list(set(genes_updated).union(set(df_merge[gene_name])))

    # cols_keep = [x for x in [gene_name, gene_ensembl_id] if x is not None]
    cols_keep = cols_gene
    cols_keep = cols_keep + cols_hgnc

    df_gene_updated = pd.concat([df[cols_keep] for df in dfs_gene.values()])
    df_m = df_gene.merge(df_gene_updated,how="left", on=cols_gene)
    df_gene_not_updated = df_m.loc[df_m["Approved symbol"].isnull()]

    if verbose:
        n_genes_found = df_gene_updated.shape[0]
        n_genes_not_found = df_gene_not_updated.shape[0]
        n_genes_updated = sum(df_gene_updated[gene_name]!=df_gene_updated["Approved symbol"])
        print("-INFO: %d gene ids were found in HGNC table" % n_genes_found)
        print("-INFO: %d gene ids were not found in HGNC table" % n_genes_not_found)
        print("-INFO: %d/%d gene symbols were updated using HGNC table" % (n_genes_updated, n_genes_found))

    df_gene = pd.concat((df_gene_updated, df_gene_not_updated))
    df_gene = df_gene.drop_duplicates()
    return df_gene


def harmonize_columns(df_gene, gene_name=None, gene_ensembl_id=None, cols_gene=[]):
    gene_name_updated = "%s_Updated" % gene_name if gene_name is not None else "Gene_Name"
    gene_ensembl_id_updated = "%s_Updated" % gene_ensembl_id if gene_ensembl_id is not None else "Ensembl_Gene_Id"
    gene_name_old = "%s_Old" % gene_name if gene_name is not None else None
    gene_ensembl_id_old = "%s_Old" % gene_ensembl_id if gene_ensembl_id is not None else None

    df_gene[gene_name_updated] = df_gene["Approved symbol"]
    df_gene[gene_ensembl_id_updated] = df_gene["Ensembl gene ID"]

    if gene_name is not None:
        mask = df_gene[gene_name_updated].isnull()
        df_gene.loc[mask, gene_name_updated] = df_gene.loc[mask, gene_name]
        df_gene = df_gene.rename(columns={gene_name: gene_name_old})
        df_gene = df_gene.rename(columns={gene_name_updated: gene_name})

    if gene_ensembl_id is not None:
        mask = df_gene[gene_ensembl_id_updated].isnull()
        df_gene.loc[mask, gene_ensembl_id_updated] = df_gene.loc[mask, gene_ensembl_id]
        df_gene = df_gene.rename(columns={gene_ensembl_id: gene_ensembl_id_old})
        df_gene = df_gene.rename(columns={gene_ensembl_id_updated: gene_ensembl_id})

    cols_keep_a = cols_gene
    cols_keep_b = [gene_name_old, gene_name, gene_ensembl_id_old, gene_ensembl_id]
    cols_keep = cols_keep_a + [x for x in cols_keep_b if x is not None and x not in cols_keep_a]
    return df_gene[cols_keep]


def main(args):
    df_gene = load_table(args.gene)
    df_hgnc = load_table(args.hgnc)

    if "." in df_gene[args.gene_ensembl_id].iloc[0]:
        df_gene[args.gene_ensembl_id + "_Version"] = df_gene[args.gene_ensembl_id]
        df_gene[args.gene_ensembl_id] = df_gene[args.gene_ensembl_id].apply(lambda x: x.split(".")[0])

    cols_gene = list(df_gene.columns)

    df_gene = update_gene_symbol_using_hgnc(df_gene=df_gene, df_hgnc=df_hgnc, gene_name=args.gene_name,
                                            gene_ensembl_id=args.gene_ensembl_id, chr_name=args.chr_name)

    df_gene = harmonize_columns(df_gene, gene_name=args.gene_name, gene_ensembl_id=args.gene_ensembl_id,
                                cols_gene=cols_gene)

    # save
    df_gene.to_csv(args.output, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update gene symbols.')
    parser.add_argument('--gene', type=str, help="Path to gene gene table.",
                        default="resources/data/gencode_v27.tsv")
    parser.add_argument('--hgnc', type=str, help='Path to hgnc table.',
                        default="../../../data/resources/hgnc/hgnc_all_symbols_03012022.tsv")
    parser.add_argument('--gene_name', type=str, default="Hugo_Symbol",
                        help="Name of the field containing gene symbols.")
    parser.add_argument('--gene_ensembl_id', type=str, default="Ensembl_Gene_Id",
                        help="Name of the field containing gene ensembl id.")
    parser.add_argument('--chr_name', type=str, default="None",
                        help="Name of the field containing the chromosome.")
    parser.add_argument('--output', type=str, help='Path to output gene table with updated symbols.')
    args = parser.parse_args()

    for arg in vars(args):
        value = getattr(args, arg)
        print("%s: %s" % (arg, value))
        if value=="None":
            setattr(args, arg, None)
    print("\n", end="")

    main(args)
