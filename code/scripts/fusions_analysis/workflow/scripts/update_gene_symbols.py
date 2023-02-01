# -*- coding: utf-8 -*-
"""
@created: 27/07/21
@modified: 27/07/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Update gene symbols in table of genes.
"""

import argparse
from functools import reduce
import numpy as np
import pandas as pd
from pyprism.data import load_table
from pyprism.util import explode_df

# function =============================================================================================================

def extract_gene_table(df, mode=None, gene_name=None, gene_ensembl_id=None, chr_name=None):
    if mode is None:
        cols = [x for x in [gene_name, gene_ensembl_id, chr_name] if x is not None]
        return df[cols].drop_duplicates()
    elif mode=="fusions":
        df_genes = []
        for i in [1,2]:
            gene_name_i = "%s_%s" % (gene_name, i)
            gene_ensembl_id_i = "%s_%s" % (gene_ensembl_id, i)
            if gene_ensembl_id_i in df:
                df_i = df[[gene_name_i, "Chr_%s" % i, gene_ensembl_id_i]].drop_duplicates()
                df_i = df_i.rename(columns={"%s_%s" % (x,i): x for x in [gene_name, "Chr", gene_ensembl_id]})
            else:
                df_i = df[[gene_name_i, "Chr_%s" % i]].drop_duplicates()
                df_i = df_i.rename(columns={"%s_%s" % (x,i): x for x in [gene_name, "Chr"]})
            df_genes.append(df_i)
        df_gene = pd.concat(df_genes).drop_duplicates()

        df_gene_cols = []
        cols_other = [x for x in ["Chr", "Gene_Id"] if "%s_1" % x in df]
        for col in cols_other:
            df_gene_col = df_gene.sort_values(by=["Gene", col])[["Gene", col]]
            df_gene_col = df_gene_col.drop_duplicates(subset="Gene", keep="first")
            agg_func = lambda x: "|".join(x.dropna().tolist())
            df_gene_col = df_gene_col.groupby("Gene").agg({col: agg_func}).reset_index()
            df_gene_col = df_gene_col.replace("", np.nan)
            df_gene_cols.append(df_gene_col)

        return reduce(lambda  left,right: pd.merge(left,right,on=['Gene'], how='outer'), df_gene_cols)


def update_gene_symbol_using_hgnc(df_gene, df_hgnc, gene_name=None, gene_ensembl_id=None, chr_name=None):
    cols_hgnc = ["Approved symbol", "Ensembl gene ID"]
    df_hgnc_approved = df_hgnc.loc[df_hgnc["Status"]=="Approved"].copy()
    df_hgnc_approved["Chr"] = df_hgnc_approved["Chromosome"].str.extract("^([0-9A-Z]+)")
    gene_ensembl_id = "Ensembl_Gene_ID_Gencode" if gene_ensembl_id is None else gene_ensembl_id

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

    cols_keep = [x for x in [gene_name, gene_ensembl_id] if x is not None]
    cols_keep = cols_keep + cols_hgnc

    df_gene_updated = pd.concat([df[cols_keep] for df in dfs_gene.values()])
    df_gene_not_updated = df_gene.loc[~df_gene[gene_name].isin(df_gene_updated[gene_name])]

    df_gene = pd.concat((df_gene_updated, df_gene_not_updated))
    df_gene = df_gene.drop_duplicates()
    return df_gene


def update_gene_ensembl_id_with_gencode(df_gene, df_gencode, gene_name=None, gene_ensembl_id=None):
    cols = df_gene.columns.tolist()

    if gene_name is not None:
        df_gencode_unique = df_gencode.sort_values(by=["tag", "gene_name"], ascending=False)
        df_gencode_unique = df_gencode_unique.drop_duplicates(subset=["gene_name"], keep="last")
        cols_gencode = ["gene_name", "Ensembl_Gene_ID_Gencode"]
        df_gene = df_gene.merge(df_gencode_unique[cols_gencode], how="left", left_on=gene_name, right_on="gene_name")

    if gene_ensembl_id is not None:
        mask = (~df_gene["Ensembl_Gene_ID_Gencode"].isnull()) & (df_gene[gene_ensembl_id].isnull())
        df_gene.loc[mask, gene_ensembl_id] = df_gene.loc[mask, "Ensembl_Gene_ID_Gencode"]
    else:
        cols += ["Ensembl_Gene_ID_Gencode"]

    return df_gene[cols]


def harmonize_columns(df_gene, gene_name=None, gene_ensembl_id=None):
    gene_name_updated = "%s_Updated" % gene_name if gene_name is not None else "Gene_Name"
    gene_ensembl_id_updated = "%s_Updated" % gene_ensembl_id if gene_ensembl_id is not None else "Ensembl_Gene_Id"

    df_gene[gene_name_updated] = df_gene["Approved symbol"]
    df_gene[gene_ensembl_id_updated] = df_gene["Ensembl gene ID"]

    if gene_name is not None:
        mask = df_gene[gene_name_updated].isnull()
        df_gene.loc[mask, gene_name_updated] = df_gene.loc[mask, gene_name]

    if gene_ensembl_id is not None:
        mask = df_gene[gene_ensembl_id_updated].isnull()
        df_gene.loc[mask, gene_ensembl_id_updated] = df_gene.loc[mask, gene_ensembl_id]

    cols_keep = [x for x in [gene_name, gene_name_updated, gene_ensembl_id, gene_ensembl_id_updated] if x is not None]
    return df_gene[cols_keep]


def update_gene_symbols(df_input, df_gene, mode, gene_name=None, gene_ensembl_id=None):
    cols = df_input.columns.tolist()
    gene_name_updated = "%s_Updated" % gene_name if gene_name is not None else "Gene_Name"
    gene_ensembl_id_updated = "%s_Updated" % gene_ensembl_id if gene_ensembl_id is not None else "Ensembl_Gene_Id"

    if mode is None:
        if gene_name is not None and gene_ensembl_id is not None:
            df_input = df_input.merge(df_gene, how="left", on=[gene_name, gene_ensembl_id])
            del df_input[gene_name]
            del df_input[gene_ensembl_id]
            df_input = df_input.rename(columns={gene_name_updated: gene_name, gene_ensembl_id_updated: gene_ensembl_id})
        elif gene_name is not None:
            df_input = df_input.merge(df_gene, how="left", on=gene_name)
            del df_input[gene_name]
            df_input = df_input.rename(columns={gene_name_updated: gene_name})
            cols += [gene_ensembl_id_updated]
        elif gene_ensembl_id is not None:
            df_input = df_input.merge(df_gene, how="left", on=gene_ensembl_id)
            del df_input[gene_ensembl_id]
            df_input = df_input.rename(columns={gene_ensembl_id_updated: gene_ensembl_id})
            cols += [gene_name_updated]

    elif mode=="fusions":
        for i in [1,2]:
            gene_name_i = "%s_%s" % (gene_name, i)
            df_input = df_input.merge(df_gene, how="left", left_on=gene_name_i, right_on=gene_name)

            del df_input[gene_name]
            del df_input[gene_name_i]

            if gene_ensembl_id is not None:
                gene_ensembl_id_i = "%s_%s" % (gene_ensembl_id, i)
                del df_input[gene_ensembl_id_i]
                gene_ensembl_id_i = "%s_%s" % (gene_ensembl_id, i)
            else:
                gene_ensembl_id_i = "%s_%s" % ("Gene_Id", i)
                cols += [gene_ensembl_id_i]

            df_input = df_input.rename(columns={gene_name_updated: gene_name_i})
            df_input = df_input.rename(columns={gene_ensembl_id_updated: gene_ensembl_id_i})

        df_input["Fusion_Id"] = df_input[["%s_1" % gene_name, "%s_2" % gene_name]].apply("--".join, axis=1)

    return df_input[cols]


def main(args):
    df_input = load_table(args.input, low_memory=False)
    df_hgnc = load_table(args.hgnc)
    df_gencode = pd.read_csv(args.gencode, sep="\t")

    df_gene = extract_gene_table(df_input, args.mode, args.gene_name, args.gene_ensembl_id, args.chr_name)
    args.chr_name = "Chr" if args.mode=="fusions" else args.chr_name
    df_gene = update_gene_ensembl_id_with_gencode(df_gene, df_gencode, args.gene_name, args.gene_ensembl_id)
    df_gene = update_gene_symbol_using_hgnc(df_gene, df_hgnc, args.gene_name, args.gene_ensembl_id, args.chr_name)
    df_gene = harmonize_columns(df_gene, args.gene_name, args.gene_ensembl_id)

    df_output = update_gene_symbols(df_input, df_gene, args.mode, args.gene_name, args.gene_ensembl_id)
    df_output.to_csv(args.output, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Update gene symbols.')
    parser.add_argument('--input', type=str, help="Path to input table with genes.",
                        default="fusions_tcga.tsv")
    parser.add_argument('--mode', type=str,
                        help="""How to extract genes from the input. Possible modes are 'fusions' or
                        None if the table is already a table of genes""",
                        default="None")
    parser.add_argument('--gene_name', type=str, default="Gene_1", help="Name of the field containing gene symbols.")
    parser.add_argument('--gene_ensembl_id', type=str, default="None",
                        help="Name of the field containing gene ensembl id.")
    parser.add_argument('--chr_name', type=str, default="Chr_1",
                        help="Name of the field containing the chromosome number.")
    parser.add_argument('--hgnc', type=str, help='Path to hgnc table.',
                        default="../../../data/resources/hgnc/hgnc_all_symbols_03012022.tsv")
    parser.add_argument('--gencode', type=str, help='Path to gencode table.',
                        default="../../../data/resources/gencode/gencode.v27.annotation_genes.tsv")
    parser.add_argument('--output', type=str, help='Path to output gene table with updated symbols.')
    args = parser.parse_args()

    if args.mode=="None":
        args.mode=None

    if args.gene_name=="None":
        args.gene_name=None

    if args.gene_ensembl_id=="None":
        args.gene_ensembl_id=None

    if args.chr_name=="None":
        args.chr_name=None

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    main(args)
