# -*- coding: utf-8 -*-
"""
@created: Dec 22 2021
@modified: Dec 22 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Take as input the tsv table with annotations produced by VEP as well as the incomplete tsv table in maf-format produced
by applying vcf2maf on the vcf file also produced by VEP.
"""

import argparse
from functools import reduce
import numpy as np
import pandas as pd

# functions ============================================================================================================

def explode_df(df, cols, sep=',', fill_value='', preserve_index=False):
    """
    Expand dataframe entries of the columns specified in l_cols and for which there are multiple values.

    Parameters
    ---------
    df: DataFrame
        Input dataframe on which expansion is performed
    cols: list
        List of columns where expansion is required
    sep: char
        Character separating the multiple values
        Default : ','
    fill_value: bool
        Entry in exanpded dataframe for empty lists
        Default : ''
    preserve_index: bool
        Whether original index should be preserved or not. If set to True, the index of the expanded DataFrame
        will be redundant.
        Default : False

    Returns
    -------
    df_expanded: DataFrame
        Returns  dataframe where entries of any of the columns in l_cols with multiple values have been expanded.
    """
    # transform comma-separated to list
    df = df.assign(**{col:df[col].str.split(sep) for col in cols}).copy()
    if (cols is not None and len(cols) > 0 and not isinstance(cols, (list, tuple, np.ndarray, pd.Series))):
        cols = [cols]
    # calculate lengths of lists
    lens = df[cols[0]].str.len()
    # format NaN to [NaN] and strip unwanted characters
    for col in cols:
        df.loc[df[col].isnull(), col] = df.loc[df[col].isnull(), col].apply(lambda x: [np.nan])
        df.loc[lens > 1, col] = df.loc[lens > 1, col].apply(lambda x: [y.strip() for y in x])
    # all columns except `cols`
    idx_cols = df.columns.difference(cols)
    # preserve original index values    
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    df_expanded = (pd.DataFrame({col:np.repeat(df[col].values, lens) for col in idx_cols},
                index=idx).assign(**{col:np.concatenate(df.loc[lens>0, col].values) for col in cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        df_expanded = (df_expanded.append(df.loc[lens==0, idx_cols], sort=False).fillna(fill_value))
    # revert the original index order
    df_expanded = df_expanded.sort_index()
    # reset index if requested
    if not preserve_index:
        df_expanded = df_expanded.reset_index(drop=True)
    return df_expanded


def identify_columns_to_be_exploded(df_vep, col_explodeby="Ensembl_transcriptid", guess_max=None):
    # Some annotations are given for all ensemble id transcripts.
    # Guess which columns are like this. Then explode them and keep only the rows where the transcript id
    # matches the annotated transcript id

    if guess_max is None:
        guess_max = df_vep.shape[0]

    df_guess = df_vep.iloc[:guess_max]
    df_guess = df_guess.replace(["-", "."], np.nan)
    df_guess = df_guess.loc[~df_guess[col_explodeby].isnull()]

    # Identify columns with commas
    cols_candidate = []
    for col in df_guess:
        vals = df_guess[col].dropna().astype(str)
        df_col_nna = df_vep.replace(["-", "."], np.nan).dropna(subset=[col])
        if any(df_col_nna[col_explodeby].isnull()):
            pass
        elif any("," in val for val in vals) and col!=col_explodeby:
            cols_candidate.append(col)

    # Identify columns that are aligned with the column col_explodeby
    cols_explode = []
    s_sizes_transcripts = df_guess[col_explodeby].apply(lambda x: len(str(x).split(",")) if x is not np.nan else 0)
    for col in cols_candidate:
        s_sizes_col = df_guess[col].apply(lambda x: len(str(x).split(",")) if x is not np.nan else 0)
        mask_nonzero = s_sizes_col!=0
        if s_sizes_transcripts[mask_nonzero].equals(s_sizes_col[mask_nonzero]):
            cols_explode.append(col)

    return cols_explode


def explode_columns(df, col_explodeby, cols_id, col_idreplace, verbose=True):
    cols_explode = identify_columns_to_be_exploded(df, col_explodeby=col_explodeby)
    cols_other = [x for x in df if x not in cols_explode]

    if verbose and len(cols_explode)>0:
        print("-the following columns will be exploded according to %s:" % col_explodeby)
        print("\t" + "\n\t".join(cols_explode))

    df_exploded = df[cols_other]
    for col in cols_explode:
        df_for_explosion = df.replace(["-", "."], np.nan)[cols_id + [col_explodeby, col]]
        df_for_explosion = df_for_explosion.dropna(subset=[col]).drop_duplicates()
        df_exploded_col = explode_df(df=df_for_explosion, cols=[col_explodeby, col], sep=",")
        df_exploded_col[col_idreplace] = df_exploded_col[col_explodeby]
        del df_exploded_col[col_explodeby]
        df_exploded_col = df_exploded_col.replace("N/A", np.nan).drop_duplicates()
        df_exploded = df_exploded.merge(df_exploded_col, how="left", on=cols_id)

    return df_exploded


def remove_empty_columns(df, name="df", ignore=[]):
    mask_empty = df.isnull().sum(axis=0)==df.shape[0]
    cols_empty = mask_empty[mask_empty].index.tolist()
    cols_empty = [x for x in cols_empty if x not in ignore]
    print("-the following empty columns are removed from %s:" % name)
    print("\t" + "\n\t".join(cols_empty))
    cols_keep = [x for x in df if x not in cols_empty]
    return df[cols_keep].copy()


def build_uploaded_variation(df_maf):
    col_id = "#Uploaded_variation"
    if df_maf.shape[0]==0:
        df_maf[col_id]=np.nan
    else:
        cols_ids = ["Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
        id_build = lambda x: "%s_%s_%s/%s" % tuple(x[cols_ids].fillna("-").values)
        df_maf[col_id] = df_maf[cols_ids].apply(id_build, axis=1)
    return df_maf.copy()


def format_value(x):
    if x is np.nan:
        return x
    else:
        is_float = type(x)==float
        x = str(x)
        if is_float and x.endswith(".0"):
            x = x.split(".0")[0]
        return x


def main(args):
    # load maf table
    df_maf = pd.read_table(args.maf_table, skiprows=1, sep="\t", na_values=["-", "."])
    del df_maf["vcf_pos"]

    # load vep table
    with open(args.vep_table, "r") as file:
        header = [x for x in file.readlines() if x.startswith("##")]
    df_vep =  pd.read_table(args.vep_table, skiprows=len(header), na_values=["-","."])
    df_vep["Feature_No_Version"] = df_vep["Feature"].apply(lambda x: x.split(".")[0] if type(x)==str else x)
    del df_vep["Allele"]

    # Some annotations are given for all ensembl id transcripts.
    # Guess which columns are like this. Then explode them and keep only the rows where the transcript id
    # matches the annotated transcript id
    col_explodeby = "Ensembl_transcriptid"
    cols_id = ["#Uploaded_variation", "Feature_No_Version"]
    col_idreplace = "Feature_No_Version"
    df_vep = explode_columns(df_vep, col_explodeby, cols_id, col_idreplace)
    del df_vep["Ensembl_transcriptid"]
    del df_vep["Ensembl_proteinid"]
    del df_vep["Feature_No_Version"]

    # remove empty columns except some
    ignore = ["#Uploaded_variation", "Chromosome", "Start_Position", "Location", "Reference_Allele",
              "Tumor_Seq_Allele2", "Feature", "HGVSc", "HGVSp", "CANONICAL"]
    df_maf = remove_empty_columns(df_maf, "df_maf", ignore=ignore)
    df_vep = remove_empty_columns(df_vep, "df_vep", ignore=ignore)

    df_vep["Chromosome"] = df_vep["Location"].apply(lambda x: x.split(":")[0])
    df_vep["Start_Position"] = df_vep["Location"].apply(lambda x: x.split(":")[1].split("-")[0])
    df_vep["Reference_Allele"] = df_vep["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[0])
    df_vep["Tumor_Seq_Allele2"] = df_vep["#Uploaded_variation"].apply(lambda x: x.split("_")[-1].split("/")[-1])
    del df_vep["Location"]
    del df_vep["#Uploaded_variation"]

    # build row identifier
    cols_ids = ["#Uploaded_variation", "Feature"]
    col_id = "Row_Identifier"
    df_maf = build_uploaded_variation(df_maf)
    df_vep = build_uploaded_variation(df_vep)
    df_maf[col_id] = df_maf[cols_ids].fillna("Unknown").apply("_".join, axis=1)
    df_vep[col_id] = df_vep[cols_ids].fillna("Unknown").apply("_".join, axis=1)
    del df_maf["#Uploaded_variation"]
    del df_vep["#Uploaded_variation"]
    df_maf = df_maf.replace(["-", "."], np.nan)
    df_vep = df_vep.replace(["-", "."], np.nan)

    # select in vep the pairs (gene, transcripts) from maf (should be CANONICAL="YES")
    df_vep = df_vep.loc[df_vep[col_id].isin(df_maf[col_id].unique())]

    # checks
    nrow_maf = df_maf.shape[0]
    nrow_vep = df_vep.shape[0]
    ids_maf = df_maf[col_id]
    ids_vep = df_vep[col_id]

    if nrow_maf!=nrow_vep:
        print("-WARNING! df_maf and df_vep have different number of rows: %d vs %d" % (nrow_maf, nrow_vep))

    if not set(df_vep["CANONICAL"].dropna().unique()).issubset(set(["YES"])):
        print("-WARNING! selected transcripts of df_vep are not all CANONICAL=YES")

    if set(ids_maf)!=set(ids_vep):
        print("-WARNING! ids_maf and ids_vep are not identical!")
        ids_maf_not_vep = set(ids_maf).difference(set(ids_vep))
        ids_vep_not_maf = set(ids_vep).difference(set(ids_maf))
        print("--ids_maf_not_vep:\n\t" + "\n\t".join(list(ids_maf_not_vep)))
        print("--ids_vep_not_maf:\n\t" + "\n\t".join(list(ids_vep_not_maf)))

    # align columns that we know have different formats
    df_vep["HGVSc"] = df_vep["HGVSc"].apply(lambda x: x.split(":")[1] if type(x)==str else x)
    df_vep["HGVSp"] = df_vep["HGVSp"].apply(lambda x: x.split(":")[1] if type(x)==str else x)

    # smart merge
    # rules:
    # 1. column names in df_maf stay
    # 2. if same colum name appear in df_maf and df_vep, assess whether they have different contents. If
    #   they do, use data from df_vep.
    # 3. for column names in df_vep and not in df_maf, assess whether this column has been renamed
    #   in df_maf. If not, add it to  df_maf.
    names_maf = list(df_maf.columns)
    names_vep = list(df_vep.columns)
    names_inter = (set(names_maf).intersection(set(names_vep))).difference(set([col_id]))
    names_maf_not_vep = set(names_maf).difference(set(names_vep))
    names_vep_not_maf = set(names_vep).difference(set(names_maf))

    # rule 2
    names_vep_replace_maf = []
    for name in names_inter:
        df_m = df_vep[[col_id, name]].merge(df_maf[[col_id, name]], how="left", on=col_id)
        name_x = "%s_x" % name
        name_y = "%s_y" % name

        # avoid mismatches caused by different types
        df_m[name_x] = df_m[name_x].apply(format_value)
        df_m[name_y] = df_m[name_y].apply(format_value)
        if not df_m[name_x].equals(df_m[name_y]):
            print("-df_maf and df_vep have different content for %s, use that of df_vep" % name)
            names_vep_replace_maf.append(name)

    # rule 3
    names_vep_add_maf = []
    for name in names_vep_not_maf:
        df_m = df_vep[[col_id, name]].merge(df_maf, how="left", on=col_id)

        match = False
        for name_maf in names_maf:
            # avoid mismatches caused by different types
            df_m[name] = df_m[name].apply(format_value)
            df_m[name_maf] = df_m[name_maf].apply(format_value)
            if df_m[name].equals(df_m[name_maf]):
                print("-df_maf and df_vep have same content for %s and %s" % (name_maf, name))
                match = True
                break

        if not match:
            names_vep_add_maf.append(name)

    # apply rules
    cols_vep = [col_id] + names_vep_replace_maf + sorted(names_vep_add_maf)
    df_maf = df_maf[[x for x in names_maf if x not in names_vep_replace_maf]]
    df_maf = df_maf.merge(df_vep[cols_vep], how="left", on=col_id)

    # delete residual columns
    del df_maf[col_id]

    # save
    df_maf.to_csv(args.output, index=False, sep="\t")

    if args.keep_vep_header:
        with open(args.output, "r") as file:
            lines = file.readlines()
        with open(args.output, "w") as file:
            for line in header:
                file.write(line)
            for line in lines:
                file.write(line)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge VEP annotations to MAF table from vcf2maf.")
    parser.add_argument('--vep_table', type=str, help='Path to tsv table produced by VEP.',
                        default="../../data/prism/wes/somatic_vep/MR424-T1-ADN_vs_MR424-N.tsv")
    parser.add_argument("--maf_table", type=str, help="Path to maf table produced by vcf2maf.",
                        default="../../data/prism/wes/somatic_vep/MR424-T1-ADN_vs_MR424-N.maf")
    parser.add_argument("--keep_vep_header", action="store_true", default=False,
                        help="If used, the header of the vep table is preserved.")
    parser.add_argument('--output', type=str, help='Path to output table.',
        default="../../data/prism/wes/somatic_maf/test.maf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
