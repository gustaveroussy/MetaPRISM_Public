# -*- coding: utf-8 -*-
"""
@created: Aug 01 2022
@modified: Oct 06 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Read current set of PASS mutations from MAF format and current filters on all mutations from VCF-like format and
apply last filtering steps.
"""

import argparse
import os
import gzip
import numpy  as     np
import pandas as     pd
import re
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


def read_table(path, **kwargs):
    header = read_header(path)
    df = pd.read_table(path, skiprows=len(header), na_values=["-","."], **kwargs)
    return df


def save_table_with_header(df, header, output):
    output_header = os.path.join(os.path.dirname(output), "header.tsv")

    print("-writing header in file %s" % output_header)
    with open(output_header, "w") as file_header:
        for line in header:
            file_header.write(line)

    # write contents
    if output.endswith(".gz"):
        output_uncompressed = output[:-3] + ".tmp"
        output_concatenate = output[:-3]
    else:
        output_uncompressed = output + ".tmp"
        output_concatenate = output

    print("-writing contents in file %s" % output_uncompressed)
    df.to_csv(output_uncompressed, index=False, sep="\t")

    # # concat both files
    cmd = "cat %s %s >> %s" % (output_header, output_uncompressed, output_concatenate)
    print("-running the command:\n\t%s" % cmd)
    subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    os.remove(output_header)
    os.remove(output_uncompressed)

    # compress if required
    if output_concatenate != output:
        cmd = "gzip %s" % output_concatenate
        print("-running the command:\n\t%s" % cmd)
        subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)



def noref_col_row_id(x, ind=-1):
    x_split = x.split("/")
    return "/".join(x_split[:ind])


def shift_col_row_id(x, shift=1, ind=-1):
    x_split = x.split("/")
    pos = x_split[ind]
    pos = int(pos) + shift
    x_split[ind] = str(pos)
    return "/".join(x_split)


def drop_duplicates_list(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def join_unique(x, sep=" || "):
    x_l = x.tolist()
    x_l = drop_duplicates_list(x_l)
    if len(x_l)>1:
        return sep.join([str(e) for e in x_l])
    else:
        return x_l[0]


def select_first_alt(df_filters):
    mask_mult_alt = df_filters["ALT"].apply(lambda x: "," in x if type(x)==str else False)
    df_filters.loc[mask_mult_alt, "ALT"] = df_filters.loc[mask_mult_alt, "ALT"].apply(lambda x: x.split(",")[0])
    return df_filters


def shift_indels_filters(df_filters):
    df_filters["Order_Shift"] = np.arange(df_filters.shape[0])

    # for simple deletions, shift POS by size of ALT
    mask_simple_del = \
        (pd.Series([x[0].startswith(x[1]) for x in zip(df_filters["REF"], df_filters["ALT"])], index=df_filters.index)) \
        & (df_filters["REF"].str.len()>df_filters["ALT"].str.len())

    # for simple insertions, shift POS by size of REF - 1
    mask_simple_ins = \
        (pd.Series([x[0].startswith(x[1]) for x in zip(df_filters["ALT"], df_filters["REF"])], index=df_filters.index)) \
        & (df_filters["ALT"].str.len()>df_filters["REF"].str.len())

    # for other indels, shift POS in df_filters by +1
    mask_indel = (df_filters["REF"].str.len()!=df_filters["ALT"].str.len())
    mask_other = mask_indel & ~mask_simple_del & ~mask_simple_ins

    # shift position
    df_filters.loc[mask_simple_del, "POS"] += df_filters.loc[mask_simple_del,"ALT"].str.len()
    df_filters.loc[mask_simple_ins, "POS"] += df_filters.loc[mask_simple_ins,"REF"].str.len()-1
    df_filters.loc[mask_other, "POS"] += 1

    # remove first nucleotide
    df_filters["REF_New"] = df_filters["REF"]
    df_filters["ALT_New"] = df_filters["ALT"]

    df_filters_a = df_filters.loc[~mask_simple_del].copy()
    df_filters_b = df_filters.loc[mask_simple_del].copy()
    remove_size = df_filters_b["ALT"].str.len()
    unique_sizes = remove_size.unique()
    for size in unique_sizes:
        mask = remove_size == size
        if sum(mask)>0:
            df_filters_b.loc[mask, "REF_New"] = df_filters_b.loc[mask, "REF"].str.slice(start=size)
            df_filters_b.loc[mask, "ALT_New"] = df_filters_b.loc[mask, "ALT"].str.slice(start=size)
    df_filters = pd.concat((df_filters_a, df_filters_b)).sort_values(by="Order_Shift")

    df_filters_a = df_filters.loc[~mask_simple_ins].copy()
    df_filters_b = df_filters.loc[mask_simple_ins].copy()
    remove_size = df_filters_b["REF"].str.len()
    unique_sizes = remove_size.unique()
    for size in unique_sizes:
        mask = remove_size == size
        if sum(mask)>0:
            df_filters_b.loc[mask, "REF_New"] = df_filters_b.loc[mask, "REF"].str.slice(start=size)
            df_filters_b.loc[mask, "ALT_New"] = df_filters_b.loc[mask, "ALT"].str.slice(start=size)
    df_filters = pd.concat((df_filters_a, df_filters_b)).sort_values(by="Order_Shift")

    df_filters.loc[mask_other, "REF_New"] = df_filters.loc[mask_other, "REF"].apply(lambda x: x[1:])
    df_filters.loc[mask_other, "ALT_New"] = df_filters.loc[mask_other, "ALT"].apply(lambda x: x[1:])
    df_filters["REF"] = df_filters["REF_New"].replace({"": np.nan})
    df_filters["ALT"] = df_filters["ALT_New"].replace({"": np.nan})
    del df_filters["REF_New"]
    del df_filters["ALT_New"]
    del df_filters["Order_Shift"]

    return df_filters


def filter_regex(df, regex, column, mode="out", return_mask=False):
    mask_regex = df[column].str.contains(regex)
    mask_regex = mask_regex.replace(np.nan, False)
    if mode=="out":
        print("-filtered out %d/%d mutations using column %s and regex %s" % \
              (sum(mask_regex), mask_regex.shape[0], column, regex))
    elif mode=="in":
        mask_regex = ~mask_regex
        print("-filtered in %d/%d mutations using column %s and regex %s" % \
              (sum(~mask_regex), mask_regex.shape[0], column, regex))
    else:
        raise ValueError("-ERROR: unrecognized value of mode %s" % mode)

    if return_mask:
        return ~mask_regex
    else:
        return df.loc[~mask_regex].copy()


def add_filter_not_exonic(df):
    col_filt = "FILTER"
    col = "Variant_Classification"
    filt = "not_exonic"
    vals_keep = ["3'UTR", "5'UTR", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Silent",
                 "Splice_Site", "Translation_Start_Site", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del",
                 "In_Frame_Ins"]
    regex_keep = "|".join(vals_keep)

    mask_keep = filter_regex(df, regex_keep, col, mode="in", return_mask=True)
    mask_pass = df[col_filt]=="PASS"
    df.loc[~mask_keep & mask_pass, col_filt] = filt
    df.loc[~mask_keep & ~mask_pass, col_filt] += ",%s" % filt

    return df


def main(args):
    # load data
    df_inp_filters = pd.read_table(args.inp_filters)
    df_inp_maf = read_table(args.inp_maf, low_memory=False)
    df_inp_maf_ann = read_table(args.inp_maf_ann, low_memory=False)
    df_inp_maf_civ = read_table(args.inp_maf_civ, low_memory=False)
    df_inp_maf_okb = read_table(args.inp_maf_okb, low_memory=False)

    # cols order
    cols_filters_order = df_inp_filters.columns.tolist()

    # add Order to preserver original order
    df_inp_filters["Order"] = np.arange(df_inp_filters.shape[0])
    df_inp_maf["Order"] = np.arange(df_inp_maf.shape[0])
    df_inp_maf_ann["Order"] = np.arange(df_inp_maf_ann.shape[0])

    # add Pair_Id
    col_pair = "Pair_Sample"
    print("-INFO: adding %s to filters, maf and maf_ann..." % col_pair)

    cols_pair = ["Tumor_Sample", "Normal_Sample"]
    df_inp_filters[col_pair] = df_inp_filters[cols_pair].fillna("NA").apply("_vs_".join, axis=1)

    cols_pair = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    df_inp_maf[col_pair] = df_inp_maf[cols_pair].fillna("NA").apply("_vs_".join, axis=1)
    df_inp_maf_ann[col_pair] = df_inp_maf_ann[cols_pair].fillna("NA").apply("_vs_".join, axis=1)

    ## 1. FILTERING on SNPs ============================================================================================

    print("-INFO: filtering on SNPs...")

    # select SNPs and create row SNP id
    col_row_id = "Row_Id"
    col_row_id_snp = "Row_Id_SNP"
    cols_filters_id = [col_pair, "CHROM", "POS", "REF", "ALT"]
    cols_maf_id = [col_pair, "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]

    print("-INFO: selecting SNPs and adding %s to filters, maf and maf_ann..." % col_row_id_snp)

    df_inp_maf_snp = df_inp_maf.loc[~df_inp_maf["Variant_Type"].isin(["INS", "DEL"])].copy()
    df_inp_maf_snp[col_row_id_snp] = df_inp_maf_snp[cols_maf_id[:-1]].astype(str).apply("/".join, axis=1)
    assert df_inp_maf_snp.shape[0]==df_inp_maf_snp[col_row_id_snp].nunique()

    print("INFO: %d SNPs in maf table" % df_inp_maf_snp.shape[0])

    df_inp_maf_ann_snp = df_inp_maf_ann.loc[~df_inp_maf_ann["Variant_Type"].isin(["INS", "DEL"])].copy()
    df_inp_maf_ann_snp[col_row_id_snp] = df_inp_maf_ann_snp[cols_maf_id[:-1]].astype(str).apply("/".join, axis=1)
    assert df_inp_maf_ann_snp.shape[0]==df_inp_maf_ann_snp[col_row_id_snp].nunique()

    print("INFO: %d SNPs in maf_ann table" % df_inp_maf_ann_snp.shape[0])

    df_inp_filters[col_row_id_snp] = df_inp_filters[cols_filters_id[:-1]].astype(str).apply("/".join, axis=1)
    df_inp_filters_snp = df_inp_filters.loc[df_inp_filters[col_row_id_snp].isin(df_inp_maf_snp[col_row_id_snp])].copy()
    df_inp_filters_oth = df_inp_filters.loc[~df_inp_filters[col_row_id_snp].isin(df_inp_maf_snp[col_row_id_snp])].copy()
    assert df_inp_filters_snp.shape[0]==df_inp_filters_snp[col_row_id_snp].nunique()

    print("INFO: %d SNPs in filters table" % df_inp_filters_snp.shape[0])

    # checks
    assert set(df_inp_maf_ann_snp[col_row_id_snp]).issubset(set(df_inp_maf_snp[col_row_id_snp]))
    assert set(df_inp_maf_snp[col_row_id_snp]).issubset(set(df_inp_filters_snp[col_row_id_snp]))

    # perform filtering
    print("INFO: adding filter not exonic on maf table...")
    df_inp_maf_snp = add_filter_not_exonic(df_inp_maf_snp)
    print("INFO: adding filter not exonic on maf_ann table...")
    df_inp_maf_ann_snp = add_filter_not_exonic(df_inp_maf_ann_snp)
    ids_pass = df_inp_maf_snp.loc[df_inp_maf_snp["FILTER"]=="PASS", col_row_id_snp].tolist()

    mask_rescue = df_inp_maf_ann_snp["FILTER"].apply(lambda x: not "not_exonic" in x)
    ids_rescue = df_inp_maf_ann_snp.loc[mask_rescue, col_row_id_snp].tolist()
    ids_pass = list(set(ids_pass).union(set(ids_rescue)))
    ids_fail = list(set(df_inp_maf_snp[col_row_id_snp]).difference(set(ids_pass)))

    # record failed ids and filtering into filters
    mask_ids_fail_filters = df_inp_filters_snp[col_row_id_snp].isin(ids_fail)
    mask_ids_fail_maf = df_inp_maf_snp[col_row_id_snp].isin(ids_fail)
    mask_ids_fail_maf_ann = df_inp_maf_ann_snp[col_row_id_snp].isin(ids_fail)
    df_inp_filters_snp_a = df_inp_filters_snp.loc[~mask_ids_fail_filters].copy()
    df_inp_filters_snp_b = df_inp_filters_snp.loc[mask_ids_fail_filters].copy()
    df_inp_maf_snp_a = df_inp_maf_snp.loc[~mask_ids_fail_maf].copy()
    df_inp_maf_ann_snp_a = df_inp_maf_ann_snp.loc[~mask_ids_fail_maf_ann].copy()
    df_inp_maf_snp_b = df_inp_maf_snp.loc[mask_ids_fail_maf].copy()

    print("INFO: merging df_inp_filters_snp_b and df_inp_maf_snp_b...")
    assert df_inp_filters_snp_b.shape[0]==df_inp_maf_snp_b.shape[0]
    del df_inp_filters_snp_b["FILTER"]
    df_inp_filters_snp_b = df_inp_filters_snp_b.merge(df_inp_maf_snp_b[[col_row_id_snp, "FILTER"]], on=col_row_id_snp)
    assert df_inp_filters_snp_b["FILTER"].isnull().sum()==0

    # reassemble with new filters
    print("INFO: concatenating df_inp_filters_snp_a and df_inp_filters_snp_b ...")
    df_inp_filters_snp = pd.concat((df_inp_filters_snp_a, df_inp_filters_snp_b))

    print("INFO: concatenating df_inp_filters_oth and df_inp_filters_snp ...")
    df_inp_filters = pd.concat((df_inp_filters_oth, df_inp_filters_snp)).sort_values(by="Order")

    ## 2. FILTERING on INDELs ==========================================================================================

    print("-INFO: filtering on INDELs...")

    # for INDELs, the coordinates and ref/alt alleles are modified when transforming the VCF file into a MAF
    # apply the changes on filters so that we can match rows afterwards
    col_row_id = "Row_Id"
    col_row_id_ind = "Row_Id_INDEL"
    cols_filters_id = [col_pair, "CHROM", "POS", "REF", "ALT"]
    cols_maf_id = [col_pair, "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]

    print("-INFO: selecting INDELs and adding %s to filters, maf and maf_ann..." % col_row_id_ind)

    # preserve original values
    df_inp_filters["POS_Original"] = df_inp_filters["POS"]
    df_inp_filters["REF_Original"] = df_inp_filters["REF"]
    df_inp_filters["ALT_Original"] = df_inp_filters["ALT"]

    # sometimes, ALT contains multiple comma-separated alternatives (when TYPE is MIXED but not only)
    # in that case, create a new "ALT" column with only 1 selected alternative
    # hypothesis: the selection rule for the maf is to keep the first alternative.
    # important: perform the selection before the shift of indels as shifts may be needed
    # note: sometimes the selection resutls in empty value for "ALT" which is an indel but for which no shift is needed
    # the function `shift_indels_filters` will not shift these.
    print("-INFO: selecting first alternative allele ...")
    df_inp_filters = select_first_alt(df_inp_filters)

    # for INDELs, the coordinates and ref/alt alleles are modified when transforming the VCF file into a MAF
    # apply the changes on filters so that we can match rows afterwards
    print("-INFO: shifting coordinates of REF/ALT ...")
    df_inp_filters = shift_indels_filters(df_inp_filters)

    df_inp_maf_ind = df_inp_maf.loc[df_inp_maf["Variant_Type"].isin(["INS", "DEL"])].copy()
    df_inp_maf_ind[col_row_id_ind] = df_inp_maf_ind[cols_maf_id].fillna("-").astype(str).apply("/".join, axis=1)
    assert df_inp_maf_ind.shape[0]==df_inp_maf_ind[col_row_id_ind].nunique()

    print("INFO: %d INDELs in maf table" % df_inp_maf_ind.shape[0])

    df_inp_maf_ann_ind = df_inp_maf_ann.loc[df_inp_maf_ann["Variant_Type"].isin(["INS", "DEL"])].copy()
    df_inp_maf_ann_ind[col_row_id_ind] = df_inp_maf_ann_ind[cols_maf_id].fillna("-").astype(str).apply("/".join, axis=1)
    assert df_inp_maf_ann_ind.shape[0]==df_inp_maf_ann_ind[col_row_id_ind].nunique()

    print("INFO: %d INDELs in maf_ann table" % df_inp_maf_ann_ind.shape[0])

    df_inp_filters[col_row_id_ind] = df_inp_filters[cols_filters_id].fillna("-").astype(str).apply("/".join, axis=1)
    df_inp_filters_ind = df_inp_filters.loc[df_inp_filters[col_row_id_ind].isin(df_inp_maf_ind[col_row_id_ind])].copy()
    df_inp_filters_oth = df_inp_filters.loc[~df_inp_filters[col_row_id_ind].isin(df_inp_maf_ind[col_row_id_ind])].copy()
    assert df_inp_filters_ind.shape[0]==df_inp_filters_ind[col_row_id_ind].nunique()

    print("INFO: %d INDELs in filters table" % df_inp_filters_ind.shape[0])

    # checks
    assert set(df_inp_maf_ann_ind[col_row_id_ind]).issubset(set(df_inp_maf_ind[col_row_id_ind]))
    assert set(df_inp_maf_ind[col_row_id_ind]).issubset(set(df_inp_filters_ind[col_row_id_ind]))

    # perform filtering
    print("INFO: adding filter not exonic on maf table...")
    df_inp_maf_ind = add_filter_not_exonic(df_inp_maf_ind)
    print("INFO: adding filter not exonic on maf_ann table...")
    df_inp_maf_ann_ind = add_filter_not_exonic(df_inp_maf_ann_ind)
    ids_pass = df_inp_maf_ind.loc[df_inp_maf_ind["FILTER"]=="PASS", col_row_id_ind].tolist()

    mask_rescue = df_inp_maf_ann_ind["FILTER"].apply(lambda x: not "not_exonic" in x)
    ids_rescue = df_inp_maf_ann_ind.loc[mask_rescue, col_row_id_ind].tolist()
    ids_pass = list(set(ids_pass).union(set(ids_rescue)))
    ids_fail = list(set(df_inp_maf_ind[col_row_id_ind]).difference(set(ids_pass)))

    # record failed ids and filtering into filters
    mask_ids_fail_filters = df_inp_filters_ind[col_row_id_ind].isin(ids_fail)
    mask_ids_fail_maf = df_inp_maf_ind[col_row_id_ind].isin(ids_fail)
    mask_ids_fail_maf_ann = df_inp_maf_ann_ind[col_row_id_ind].isin(ids_fail)
    df_inp_filters_ind_a = df_inp_filters_ind.loc[~mask_ids_fail_filters].copy()
    df_inp_filters_ind_b = df_inp_filters_ind.loc[mask_ids_fail_filters].copy()
    df_inp_maf_ind_a = df_inp_maf_ind.loc[~mask_ids_fail_maf].copy()
    df_inp_maf_ann_ind_a = df_inp_maf_ann_ind.loc[~mask_ids_fail_maf_ann].copy()
    df_inp_maf_ind_b = df_inp_maf_ind.loc[mask_ids_fail_maf].copy()

    print("INFO: merging df_inp_filters_ind_b and df_inp_maf_ind_b...")
    assert df_inp_filters_ind_b.shape[0]==df_inp_maf_ind_b.shape[0]
    del df_inp_filters_ind_b["FILTER"]
    df_inp_filters_ind_b = df_inp_filters_ind_b.merge(df_inp_maf_ind_b[[col_row_id_ind, "FILTER"]], on=col_row_id_ind)
    assert df_inp_filters_ind_b["FILTER"].isnull().sum()==0

    # reassemble with new filters
    print("INFO: concatenating df_inp_filters_ind_a and df_inp_filters_ind_b ...")
    df_inp_filters_ind = pd.concat((df_inp_filters_ind_a, df_inp_filters_ind_b))

    print("INFO: concatenating df_inp_filters_oth and df_inp_filters_ind ...")
    df_inp_filters = pd.concat((df_inp_filters_oth, df_inp_filters_ind)).sort_values(by="Order")

    # SNP and INDEL assembly and remove added cols =====================================================================

    print("INFO: concatenating df_inp_maf_snp_a and df_inp_maf_snp_b ...")
    df_inp_maf = pd.concat((df_inp_maf_snp_a, df_inp_maf_ind_a), axis=0).sort_values(by="Order")

    print("INFO: concatenating df_inp_maf_ann_snp_a and df_inp_maf_ann_snp_b ...")
    df_inp_maf_ann = pd.concat((df_inp_maf_ann_snp_a, df_inp_maf_ann_ind_a), axis=0).sort_values(by="Order")
    df_inp_filters = df_inp_filters.sort_values(by="Order")

    del df_inp_maf["Order"]
    del df_inp_maf[col_pair]
    del df_inp_maf[col_row_id_snp]
    del df_inp_maf[col_row_id_ind]

    del df_inp_maf_ann["Order"]
    del df_inp_maf_ann[col_pair]
    del df_inp_maf_ann[col_row_id_snp]
    del df_inp_maf_ann[col_row_id_ind]

    del df_inp_filters["Order"]
    del df_inp_filters[col_pair]
    del df_inp_filters[col_row_id_snp]
    del df_inp_filters[col_row_id_ind]
    del df_inp_filters["POS"]
    del df_inp_filters["REF"]
    del df_inp_filters["ALT"]
    df_inp_filters = df_inp_filters.rename(columns={"POS_Original": "POS", "REF_Original": "REF",
                                                    "ALT_Original": "ALT"})

    df_inp_filters = df_inp_filters[cols_filters_order]

    # save refiltered tables
    print("INFO: saving df_inp_filters to %s ..." % args.out_filters)
    df_inp_filters.to_csv(args.out_filters, sep="\t", index=False)

    print("INFO: saving df_inp_maf to %s ..." % args.out_maf)
    header_inp_maf = read_header(args.inp_maf)
    save_table_with_header(df_inp_maf, header_inp_maf, args.out_maf)

    print("INFO: saving df_inp_maf_ann to %s ..." % args.out_maf_ann)
    df_inp_maf_ann.to_csv(args.out_maf_ann, sep="\t", index=False)

    # for CIV and OKB, only remove non-exonic ==========================================================================
    print("INFO: adding filter not exonic on maf civ table...")
    df_inp_maf_civ = add_filter_not_exonic(df_inp_maf_civ)
    mask_rescue = df_inp_maf_civ["FILTER"].apply(lambda x: not "not_exonic" in x)
    df_inp_maf_civ = df_inp_maf_civ.loc[mask_rescue].copy()
    print("INFO: removed %d/%d non-exonic variants from maf civic table" % (sum(~mask_rescue), len(mask_rescue)))

    print("INFO: adding filter not exonic on maf okb table...")
    df_inp_maf_okb = add_filter_not_exonic(df_inp_maf_okb)
    mask_rescue = df_inp_maf_okb["FILTER"].apply(lambda x: not "not_exonic" in x)
    df_inp_maf_okb = df_inp_maf_okb.loc[mask_rescue].copy()
    print("INFO: removed %d/%d non-exonic variants from maf oncokb table" % (sum(~mask_rescue), len(mask_rescue)))

    print("INFO: saving df_inp_maf_civ to %s ..." % args.out_maf_civ)
    df_inp_maf_civ.to_csv(args.out_maf_civ, sep="\t", index=False)

    print("INFO: saving df_inp_maf_okb to %s ..." % args.out_maf_okb)
    df_inp_maf_okb.to_csv(args.out_maf_okb, sep="\t", index=False)


if __name__ == "__main__":
    # parameters =======================================================================================================

    parser = argparse.ArgumentParser(description='Apply last mutation filters.')
    parser.add_argument('--inp_filters', type=str, help='Path to all mutations with current filtering.')
    parser.add_argument('--inp_maf', type=str, help='Path to all currently PASS mutations.')
    parser.add_argument('--inp_maf_ann', type=str, help='Path to all currently PASS and oncogenic mutations.')
    parser.add_argument('--inp_maf_civ', type=str, help='Path to all currently PASS and CIViC-annotated mutations.')
    parser.add_argument('--inp_maf_okb', type=str, help='Path to all currently PASS and OncoKB-annotated mutations.')
    parser.add_argument('--out_filters', type=str, help='Path to all mutations with new filtering.')
    parser.add_argument('--out_maf', type=str, help='Path to all newly PASS mutations.')
    parser.add_argument('--out_maf_ann', type=str, help='Path to all newly PASS and oncogenic mutations.')
    parser.add_argument('--out_maf_civ', type=str, help='Path to all newly PASS and CIViC-annotated mutations.')
    parser.add_argument('--out_maf_okb', type=str, help='Path to all newly PASS and OncoKB-annotated mutations.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")
    main(args)
