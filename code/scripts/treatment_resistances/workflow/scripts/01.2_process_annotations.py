# -*- coding: utf-8 -*-
"""
@created: Aug 08 2021
@modified: Dec 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Convert manual annotations of treatment resistances to a computer-friendly format ready for analysis.
"""

import argparse
from functools import reduce
import numpy as np
import pandas as pd
import re
import sys

# pyprism
from pyprism.data import load_table, load_cln, split_targeted_therapy_targets
from pyprism.util import explode_df

def select_subjects_and_samples(df, df_cln, col_tt, col_sub_id="Subject_Id", tumor_type="All",
                                sample_type="DNA_T\\|RNA_T", dnap=True):
    # select patients with from tumor type
    if tumor_type == "All":
        subjects = df_cln[col_sub_id].unique().tolist()
    else:
        mask_keep = df_cln[col_tt].isin(tumor_type.split("__"))
        subjects = df_cln.loc[mask_keep,col_sub_id].unique().tolist()

    mask = df[col_sub_id].isin(subjects)
    df_sub = df.loc[mask].copy()
    print("-selected %d/%d lines from df from selection of %d subjects" % (sum(mask), len(mask), len(subjects)))

    # select patients with selected data
    mask_keep = df_cln["Sample_Type"].apply(lambda x: bool(re.search(sample_type, x)))
    subjects = df_cln.loc[mask_keep, col_sub_id].unique()
    if dnap:
        samples_dna = df_cln.loc[mask_keep, "Sample_Id_DNA_P"].unique().tolist()
    else:
        samples_dna = df_cln.loc[mask_keep, "Sample_Id_DNA_T"].unique().tolist()
    samples_rna = df_cln.loc[mask_keep, "Sample_Id_RNA_T"].unique().tolist()

    mask = df_sub[col_sub_id].isin(subjects)
    df_sub = df_sub.loc[mask].copy()
    print("-selected %d/%d lines from df from selection of %d subjects" % (sum(mask), len(mask), len(subjects)))

    return df_sub


def compare_identifiers(df_a, df_b, col_id, cols_info=[], name_a="A", name_b="B"):
    ids_a = set(df_a[col_id])
    ids_b = set(df_b[col_id])
    ids_a_not_b = ids_a.difference(ids_b)
    ids_b_not_a = ids_b.difference(ids_a)
    df_a_not_b = df_a.loc[df_a[col_id].isin(list(ids_a_not_b))]
    df_b_not_a = df_b.loc[df_b[col_id].isin(list(ids_b_not_a))]

    print("-INFO: %s ids in %s and not in %s" % (len(ids_a_not_b), name_a, name_b))
    print("-INFO: %s ids in %s and not in %s" % (len(ids_b_not_a), name_b, name_a))

    for col_info in cols_info:
        print("-here are the counts of %s for ids in %s not in %s" % (col_info, name_a, name_b))
        print(df_a_not_b[col_info].value_counts())

        print("-here are the counts of %s for ids in %s not in %s" % (col_info, name_b, name_a))
        print(df_b_not_a[col_info].value_counts())


def standardize_column(x, sep=","):
    x_split = x.split(sep)
    x_split = [e.strip().upper() for e in x_split]
    return sep.join(x_split)


def clean_reference(x):
    if type(x)==float:
        return x
    else:
        x = x.replace(u'\xa0', u' ')
        x = re.sub("^pmid ", "PMID: ", x)
        x = re.sub("^PMC", "PMC: ", x)
        x = re.sub("^PMC: :", "PMC: ", x)
        x = re.sub("^doi", "DOI", x)
    return x


def unique_list_keep_order(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def join_drop_na(x, sep="|"):
    x_nna = x.replace({"N/A": np.nan}).dropna().tolist()
    x_nna = unique_list_keep_order(x_nna)
    if len(x_nna)==0:
        return "N/A"
    else:
        return sep.join(x_nna)

def main(args):
    # variables
    cohort = "prism"
    col_row_id = "Row_Identifier"
    col_gen = "Hugo_Symbol"
    col_alt = "Alteration"
    col_alt_cat = "Alteration_Category"
    col_alt_det = "Alteration_Detail"
    col_sub_id = "Subject_Id"
    col_tt = "Tumor_Type"
    col_tt_okb = "MSKCC_Oncotree"
    col_tt_civ = "Civic_Disease"
    col_drug = "Treatments"

    # load data
    df_cln = load_cln(cohort)
    df_cln = df_cln.rename(columns={"Project_TCGA_More": col_tt, "Drugs_Before_Biopsy": col_drug})
    df_cln["Sample_Id_DNA_P"] = df_cln[["Sample_Id_DNA_T", "Sample_Id_DNA_N"]].fillna("NA").apply("_vs_".join, axis=1)
    df_alt_bef = load_table(args.alterations_before_annot)
    df_alt_aft = load_table(args.alterations_after_annot)
    df_alt_cur = load_table(args.alterations_current)

    # for DNA, replace DNA_T by DNA_P
    df_ids = df_cln.loc[:,["Sample_Id_DNA_T", "Sample_Id_DNA_P"]].dropna()
    df_ids = df_ids.rename(columns={"Sample_Id_DNA_T": "Sample_Id"})
    df_alt_bef = df_alt_bef.merge(df_ids, how="left", on="Sample_Id")
    df_alt_bef["Sample_Id"]
    df_alt_bef_dna = df_alt_bef.loc[~df_alt_bef["Alteration_Category"].isin(["Fusion"])].copy()
    df_alt_bef_dna["Sample_Id_DNA_P"].isnull().sum()
    df_alt_bef_dna.loc[df_alt_bef_dna["Sample_Id_DNA_P"].isnull(), "Sample_Id"].unique()

    # manual corrections
    df_alt_aft["Resistances"] = df_alt_aft["Resistances"].replace({0: np.nan})
    df_alt_aft["Reference"]  = df_alt_aft["Reference"].apply(clean_reference)

    # apply relabelling
    for df in [df_alt_bef, df_alt_aft]:
        for old, new in zip(["TMB-High", "MSI-High"], ["TMB High", "MSI High"]):
            df[col_row_id] = df[col_row_id].apply(lambda x: x.replace(old, new))
            df[col_alt] = df[col_alt].replace({old: new})

    # select samples
    tumor_type = "__".join(df_alt_bef[col_tt].unique().tolist())
    df_alt_bef = select_subjects_and_samples(df_alt_bef, df_cln, col_tt, col_sub_id=col_sub_id, tumor_type=tumor_type,
                                             sample_type="DNA_T\\|RNA_T", dnap=False)
    df_alt_aft = select_subjects_and_samples(df_alt_aft, df_cln, col_tt, col_sub_id=col_sub_id, tumor_type=tumor_type,
                                             sample_type="DNA_T\\|RNA_T", dnap=False)
    df_alt_cur = select_subjects_and_samples(df_alt_cur, df_cln, col_tt, col_sub_id=col_sub_id, tumor_type=tumor_type,
                                             sample_type="DNA_T\\|RNA_T")

    # compare tables
    # compare_identifiers(df_a=df_alt_bef, df_b=df_alt_cur, col_id=col_row_id, name_a="Before", name_b="Current",
    #                     cols_info=[col_alt_cat, "Variant_Classification", col_tt])

    assert set(df_alt_bef[col_row_id]).issubset(set(df_alt_bef[col_row_id]))

    # compare_identifiers(df_a=df_alt_bef, df_b=df_alt_aft, col_id=col_row_id, name_a="Before", name_b="After",
    #                     cols_info=[col_alt_cat, col_tt])

    assert set(df_alt_aft[col_row_id]).issubset(set(df_alt_bef[col_row_id]))

    # for alterations after annotations, keep only annotated lines
    cols_ann = [x for x in  df_alt_aft if x not in df_alt_bef]
    df_alt_aft = df_alt_aft.loc[~df_alt_aft["Resistances"].isnull()].copy()
    df_alt_aft["Resistances"] = df_alt_aft["Resistances"].apply(standardize_column)

    # format "Resistances" entries to have standard values
    ids_resist_a = set(df_alt_aft[col_row_id])
    df_alt_aft = explode_df(df=df_alt_aft, cols=["Resistances"], sep=",")
    df_alt_aft = explode_df(df=df_alt_aft, cols=["Treatments"], sep="|")
    df_alt_aft = df_alt_aft.loc[df_alt_aft["Treatments"]==df_alt_aft["Resistances"]]
    ids_resist_b = set(df_alt_aft[col_row_id])
    assert ids_resist_b==ids_resist_a

    cols_keep = [col_tt, col_tt_okb, col_tt_civ, col_gen, col_alt, col_alt_cat] \
            + ["Resistances", "Reference"]
    df_alt_aft_rules = df_alt_aft[cols_keep].fillna("N/A").groupby(cols_keep).size().to_frame("Count").reset_index()
    df_alt_aft_rules = df_alt_aft_rules.sort_values(by=cols_keep)

    # integrate manual review of rules
    cols_merge_id_rev = [col_tt, col_tt_okb, col_tt_civ, col_gen, col_alt, "Resistances"]
    col_merge_id_rev = "Merge_Id_Rev"
    filepath_rev = "../../../results/treatment_resistances/annotations/rules_annotation_LF_review_SN.xlsx"

    df_alt_aft_rules_rev = pd.read_excel(filepath_rev)
    df_alt_aft_rules[col_merge_id_rev] = df_alt_aft_rules[cols_merge_id_rev].apply("/".join, axis=1)
    df_alt_aft_rules_rev[col_merge_id_rev] = df_alt_aft_rules_rev[cols_merge_id_rev].apply("/".join, axis=1)

    mask_a = df_alt_aft_rules[col_merge_id_rev].isin(df_alt_aft_rules_rev[col_merge_id_rev])
    df_alt_aft_rules_a = df_alt_aft_rules.loc[~mask_a].copy()
    df_alt_aft_rules_b = df_alt_aft_rules.loc[mask_a].copy()

    df_alt_aft_rules_rev = df_alt_aft_rules_rev.loc[df_alt_aft_rules_rev["Reference"]!="exclude"].copy()
    df_alt_aft_rules = pd.concat((df_alt_aft_rules_a, df_alt_aft_rules_rev))
    del df_alt_aft_rules[col_merge_id_rev]
    df_alt_aft_rules.to_excel(args.output_rules_annot, index=False)

    # apply rules to all alterations
    cols_merge_id = [col_tt, col_tt_okb, col_tt_civ, col_gen, col_alt]
    col_merge_id = "Merge_Id_Annot"
    df_alt_aft_rules[col_merge_id] = df_alt_aft_rules[cols_merge_id].apply("/".join, axis=1)

    df_alt_aft_rules = df_alt_aft_rules.groupby(col_merge_id).agg({"Resistances": lambda x: join_drop_na(x, "|"),
                                                                   "Reference": lambda x: join_drop_na(x, ";")})
    df_alt_aft_rules = df_alt_aft_rules.reset_index()
    df_alt_aft_rules = df_alt_aft_rules.rename(columns={"Resistances": "LF_Res_Drug", "Reference": "LF_Res_Citations"})

    df_alt_cur[col_merge_id] = df_alt_cur[cols_merge_id].apply("/".join, axis=1)
    df_alt_cur = df_alt_cur.merge(df_alt_aft_rules, how="left", on=col_merge_id)
    del df_alt_cur[col_merge_id]
    df_alt_cur.to_csv(args.output_alterations, sep="\t", index=False)
    print("-table saved at %s" % args.output_alterations)

    # make table with 1 line per drug received and columns describing mechanisms of resistance
    df_cln = select_subjects_and_samples(df_cln, df_cln, col_tt, col_sub_id=col_sub_id, tumor_type=tumor_type,
                                         sample_type="DNA_T\\|RNA_T")
    cols_keep = [col_sub_id, col_tt, col_tt_okb, col_tt_civ, col_drug]
    df_cln = df_cln[cols_keep].copy()
    mask_drop = (~df_cln[col_drug].isnull()) & (df_cln[col_drug]!="None")
    print("-selecting %d/%d subjects from df_cln with no or missing treatment" % (sum(mask_drop), len(mask_drop)))
    df_cln = df_cln.loc[mask_drop].copy()

    # prepare 1 line per drug
    df_cln = explode_df(df_cln, cols=[col_drug], sep="|")

    # add columns for OncoKB/CIViC/LF
    col_evt = "Gene_Alteration"
    df_alt_cur[col_evt] = df_alt_cur[[col_gen, col_alt]].apply(" ".join, axis=1)

    dfs_cln_tag = []
    for tag in ["Oncokb", "Civic", "LF"]:
        col_lvl = "%s_Res_Level_Simple" % tag
        col_res = "%s_Res_Drug" % tag
        col_cite = "%s_Res_Citations" % tag
        cols_tag = [x for x in [col_sub_id, col_evt, col_lvl, col_res, col_cite] if x in df_alt_cur]
        df_alt_cur_tag = df_alt_cur[cols_tag].copy()
        df_alt_cur_tag = df_alt_cur_tag.loc[~df_alt_cur_tag[col_res].isnull()].copy()
        df_alt_cur_tag = df_alt_cur_tag.rename(columns={col_res: col_drug})
        df_alt_cur_tag = explode_df(df=df_alt_cur_tag, cols=[col_drug], sep="|")
        df_alt_cur_tag = df_alt_cur_tag.drop_duplicates()
        df_cln_tag = df_cln.merge(df_alt_cur_tag, how="left", on=[col_sub_id, col_drug])
        cols_agg = [x for x in [col_evt, col_lvl, col_cite] if x in df_alt_cur_tag]
        df_cln_tag = df_cln_tag.fillna("N/A")
        df_cln_tag = df_cln_tag.groupby([col_sub_id, col_drug]).agg({col_agg: lambda x: join_drop_na(x, " & ") \
                                                                     for col_agg in cols_agg})
        df_cln_tag = df_cln_tag.reset_index(drop=False)
        df_cln_tag = df_cln_tag.rename(columns={col_evt: "%s_%s" % (tag, col_evt)})
        dfs_cln_tag.append(df_cln_tag)

    df_cln_trt = reduce(lambda df1,df2: pd.merge(df1,df2,on=[col_sub_id, col_drug]), [df_cln] + dfs_cln_tag)
    df_cln_trt = df_cln_trt.replace({"N/A": np.nan})
    df_cln_trt.to_excel(args.output_treatments, index=False)
    print("-table saved at %s" % args.output_treatments)


# run ==================================================================================================================

if __name__ == "__main__":
    results_folder = "../../../results/treatment_resistances"

    parser = argparse.ArgumentParser(description="Oncoplot-like figure detailing treatment resistances.")
    parser.add_argument('--alterations_before_annot', type=str, help='Path to table of alterations before annotations.',
                        default="%s/annotations/aggregated_alterations_prism_072022.xlsx" % results_folder)
    parser.add_argument('--alterations_after_annot', type=str, help='Path to table of alterations after annotations.',
                        default="%s/annotations/aggregated_alterations_prism_annotated_LF_072022.xlsx" % results_folder)
    parser.add_argument('--alterations_current', type=str, help='Path to table of current alterations.',
                        default="../../../results/combined_alterations/alterations/aggregated_alterations_prism_all.tsv")
    parser.add_argument('--output_rules_annot', type=str, help='Path to table annotation rules.',
                        default="%s/annotations/rules_annotation_LF.xlsx" % results_folder)
    parser.add_argument('--output_alterations', type=str,
                        help='Path to table of alterations with additonal annotations.',
                        default="%s/annotations/aggregated_alterations_prism_all_more.tsv" % results_folder)
    parser.add_argument('--output_treatments', type=str,  help='Path to table with 1 line per treatment received.',
                        default="%s/annotations/treatment_resistances_prism.xlsx" % results_folder)
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
