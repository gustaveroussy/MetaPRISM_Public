# -*- coding: utf-8 -*-
"""
@created: 01 Feb 22
@modified: 02 Jan 23
@authors: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Aggregate mutations annotations from CIViC, OncoKB and CGI (to be done).
"""

import argparse
import gzip
import numpy as np
import pandas as pd
import re
import sys

# functions ============================================================================================================

def merge_on_rows(df_x, df_y, how_rows="inner", prefer_y_over_x=True, prefer_exceptions: list=[],
                  fillna_where_possible=False, filepath_discrepancies: str=None, name_x: str="X", name_y: str="Y",
                  remove_na_from_discrepancies=False, write_only_discrepant_cols=False):
    rows_x = set(df_x.index)
    rows_y = set(df_y.index)
    rows_common = rows_x.intersection(rows_y)
    rows_x_only = rows_x.difference(set(rows_common))
    rows_y_only = rows_y.difference(set(rows_common))
    df_xn = df_x.copy()
    df_yn = df_y.copy()
    n_common = len(rows_common)

    rows_common = list(rows_common)
    rows_x_only = list(rows_x_only)
    rows_y_only = list(rows_y_only)

    if list(df_xn.index.names) == [None] or list(df_yn.index.names) == [None]:
        df_xn.index.names = ["Internal_Index"]
        df_yn.index.names = ["Internal_Index"]

    cols_x = set(df_xn.columns)
    cols_y = set(df_yn.columns)
    cols_common = cols_x.intersection(cols_y)
    cols_x_only = cols_x.difference(set(cols_common))
    cols_y_only = cols_y.difference(set(cols_common))

    def _add_level_columns(df, name, value, pos=0):
        df_cols = df.columns.to_frame()
        df_cols.insert(pos, name, value)
        df.columns = pd.MultiIndex.from_frame(df_cols)
        return df

    cols_common_drop = []
    dfs_discrepant = {}
    for col_c in cols_common:
        df_x_c = df_xn.loc[rows_common, col_c].fillna("NaN")
        df_y_c = df_yn.loc[rows_common, col_c].fillna("NaN")
        mask_discrepant = df_x_c != df_y_c
        n_discrepant = sum(mask_discrepant)

        if n_discrepant != 0:
            if remove_na_from_discrepancies:
                mask_discrepant = mask_discrepant & ~((df_x_c=="NaN") | (df_y_c=="NaN"))

            if sum(mask_discrepant) > 0:
                print("-warning: %s is discrepant for %d/%d common individuals" % \
                      (col_c, sum(mask_discrepant), n_common), flush=True)
            cols_common_drop.append(col_c)

            df_discrepant_x = df_xn.loc[rows_common].loc[mask_discrepant]
            if write_only_discrepant_cols:
                df_discrepant_x = df_discrepant_x[[col_c]]
            else:
                df_discrepant_x = df_discrepant_x[[c for c in df_discrepant_x.columns if c!=col_c] + [col_c]]
            df_discrepant_x = _add_level_columns(df_discrepant_x, name_x, name_x)

            df_discrepant_y = df_yn.loc[rows_common].loc[mask_discrepant]
            if write_only_discrepant_cols:
                df_discrepant_y = df_discrepant_y[[col_c]]
            else:
                df_discrepant_y = df_discrepant_y[[col_c] + [c for c in df_discrepant_y.columns if c!=col_c]]
            df_discrepant_y = _add_level_columns(df_discrepant_y, name_y, name_y)
            df_discrepant = pd.concat((df_discrepant_x, df_discrepant_y), axis=1).sort_index()

            dfs_discrepant["%s" % str(col_c)] = df_discrepant

    if filepath_discrepancies is not None and any([df.shape[0] > 0 for df in dfs_discrepant.values()]):
        with pd.ExcelWriter(filepath_discrepancies) as writer:
            for name, df in sorted(dfs_discrepant.items()):
                if df.shape[0] > 0:
                    one_saved = True
                    save_df_to_data(df, writer, sheet_name=name[:31], index=True, header=True, verbose=False)
            print("-excel file with discrepancies saved at %s" % filepath_discrepancies, flush=True)

    for col_c in cols_common_drop:
        if  col_c in prefer_exceptions:
            choose_y = not prefer_y_over_x
        else:
            choose_y = prefer_y_over_x

        if choose_y:
            del df_xn[col_c]
        else:
            del df_yn[col_c]
        cols_common.remove(col_c)

    index_names = df_xn.index.names
    df_xn = df_xn.reset_index(drop=False)
    df_yn = df_yn.reset_index(drop=False)
    cols_common = list(cols_common) + list(index_names)
    df_xn.loc[:, cols_common] = df_xn[cols_common].fillna("NaN")
    df_yn.loc[:, cols_common] = df_yn[cols_common].fillna("NaN")
    df_c = df_xn.merge(df_yn, how=how_rows, on=cols_common)
    df_c = df_c.replace(to_replace="NaN", value=np.nan)
    df_c = df_c.set_index(index_names)

    if list(df_c.index.names) == ["Internal_Index"]:
        df_c.index.names = [None]

    # fill in the gaps where possible
    if fillna_where_possible:
        for col_c in cols_common_drop:
            rows_na_x = df_x.loc[df_x[col_c].isnull()].index
            rows_common_na_x = list(set(rows_common).intersection(set(rows_na_x)))
            df_c.loc[rows_common_na_x, col_c] = df_y.loc[rows_common_na_x, col_c]

            rows_na_y = df_y.loc[df_y[col_c].isnull()].index
            rows_common_na_y = list(set(rows_common).intersection(set(rows_na_y)))
            df_c.loc[rows_common_na_y, col_c] = df_x.loc[rows_common_na_y, col_c]

    if len(rows_x_only) > 0:
        for col_x in cols_x:
            df_c.loc[rows_x_only, col_x]  = df_x.loc[rows_x_only, col_x]

    if len(rows_y_only) > 0:
        for col_y in cols_y:
            df_c.loc[rows_y_only, col_y]  = df_y.loc[rows_y_only, col_y]

    return df_c


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


def build_row_identifier_maf(df, col_gene="Hugo_Symbol", col_start="Start_Position", col_end="End_Position",
                             col_ref="Reference_Allele", col_alt="Tumor_Seq_Allele2"):
    dfc = df.copy()
    cols = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", col_gene, col_start, col_end, col_ref, col_alt]
    dfc[col_start] =  dfc[col_start].apply(lambda x: "%d" % x if not np.isnan(x) else "")
    dfc[col_end] =  dfc[col_end].apply(lambda x: "%d" % x if not np.isnan(x) else "")
    cols = [x for x in cols if x in dfc]
    dfc["Row_Identifier"] = dfc[cols].fillna("-").astype(str).apply("_".join, axis=1)
    df["Row_Identifier"] = dfc["Row_Identifier"]
    return df


def build_row_identifier_cna(df, col_gene="Hugo_Symbol", col_alt="Alteration"):
    dfc = df.copy()
    cols = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", col_gene, col_alt]
    cols = [x for x in cols if x in dfc]
    dfc["Row_Identifier"] = dfc[cols].fillna("-").astype(str).apply("_".join, axis=1)
    df["Row_Identifier"] = dfc["Row_Identifier"]
    return df


def build_row_identifier(df, cat, **kwargs):
    if args.cat=="maf":
        return build_row_identifier_maf(df, **kwargs)
    elif args.cat=="cna":
        return build_row_identifier_cna(df, **kwargs)


def correct_dtype_from_num_to_str(x):
    try:
        y = "%d" % int(float(x))
    except:
        try:
            y = str(x)
            if y=="nan":
                y = np.nan
        except:
            y = x

    return y


def add_pair_id(df):
    cols_a = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    cols_b = ["Sample_Id_DNA_T", "Sample_Id_DNA_N"]
    if set(cols_a).issubset(set(df.columns)):
        cols = cols_a
    elif set(cols_b).issubset(set(df.columns)):
        cols = cols_b
    else:
        raise ValueError("-no column available for getting pair id")
    df["Pair_Id"] = df[cols].fillna("NA").apply("_".join, axis=1)
    return df


def main(args):
    df_okb = pd.read_table(args.okb)
    df_civ = pd.read_table(args.civ)

    # reorder columns civ
    cols_ppd = [x for x in df_civ if re.search("^Prognostic:|^Predictive:|^Diagnostic:", x) is not None]
    cols_oth = [x for x in df_civ if x not in cols_ppd]
    df_civ = df_civ[cols_oth+cols_ppd]

    if args.cat == "cna":
        df_okb = df_okb.rename(columns={"HUGO_SYMBOL": "Hugo_Symbol", "ALTERATION": "Alteration"})

    df_okb = build_row_identifier(df_okb, args.cat)
    df_civ = build_row_identifier(df_civ, args.cat)

    if args.cat == "maf":
        # harmonize common columns
        df_okb["NCBI_Build"] = "hg19"

    for df in [df_okb, df_civ]:
        cols = df.columns.tolist()
        for col in cols:
            df[col] = df[col].apply(correct_dtype_from_num_to_str)

    # info messages
    okb_and_civ = len(set(df_okb["Row_Identifier"]).intersection(set(df_civ["Row_Identifier"])))
    okb_not_civ = len(set(df_okb["Row_Identifier"]).difference(set(df_civ["Row_Identifier"])))
    civ_not_okb = len(set(df_civ["Row_Identifier"]).difference(set(df_okb["Row_Identifier"])))
    okb_or_civ = len(set(df_civ["Row_Identifier"]).union(set(df_okb["Row_Identifier"])))
    print("-INFO: there are %d rows in OncoKB and CIViC" % okb_and_civ)
    print("-INFO: there are %d rows in OncoKB not in CIViC" % okb_not_civ)
    print("-INFO: there are %d rows in CIViC not in OncoKB" % civ_not_okb)
    print("-INFO: there are %d rows in OncoKB or CIViC" % okb_or_civ)

    df_ann = merge_on_rows(df_x=df_okb.set_index("Row_Identifier"),
                           df_y=df_civ.set_index("Row_Identifier"),
                           how_rows="outer",
                           prefer_y_over_x=True,
                           prefer_exceptions=[],
                           fillna_where_possible=True,
                           filepath_discrepancies=None,
                           name_x="OncoKB",
                           name_y="CIViC",
                           remove_na_from_discrepancies=True,
                           write_only_discrepant_cols=True)

    df_ann = df_ann.reset_index(drop=True)

    # check
    assert df_ann.shape[0] == okb_or_civ
    df_ann = df_ann.replace("nan", np.nan)

    # for mutations, Civic annotations are sloppy. Retain only annotations for mutations also 
    # annotated in Oncokb
    if args.cat == "maf":
        if "ONCOGENIC" in df_ann:
            mask_okb = ~df_ann["ONCOGENIC"].isnull()
        else:
            mask_okb = pd.Series(False, index=df_ann.index)

        if "CIViC_Matching_Gene_Variant" in df_ann:
            mask_civ = ~df_ann["CIViC_Matching_Gene_Variant"].isnull()
        else:
            mask_civ = pd.Series(False, index=df_ann.index)

        mask_civ_not_okb = mask_civ & ~mask_okb
        df_ann = df_ann.loc[~mask_civ_not_okb].copy()

    # save
    df_ann.to_csv(args.output, index=False, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate mutations annotations.")
    parser.add_argument("--civ", type=str, help="Path to CIViC-annotated mutations table.")
    parser.add_argument("--okb", type=str, help="Path to OncoKB-annotated mutations table.")
    parser.add_argument("--cat", type=str, help="Category of annotations. One of 'maf' or 'cna'.", default='maf')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
