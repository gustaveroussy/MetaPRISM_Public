# -*- coding: utf-8 -*-
"""
@created: Feb 24 2022
@modified: Feb 28 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France
"""

import argparse
import os
import numpy as np
import pandas as pd
import re
import sys
from tqdm import tqdm
sys.path.append("workflow/functions")

import seaborn as sns
import matplotlib.pyplot as plt

# set seaborn whitegrid theme
sns.set(style="whitegrid")

from filter import filter_out_fusions, filter_in_whitelist
from utils import load_table_samples, add_aggregated_callers
from pyprism.data import load_bio, load_table

# functions ============================================================================================================

def load_table_samples(algo, cohort):
    filepath = "../../../data/%s/rna/%s/sample_list.tsv" % (cohort, algo)
    df_sam = pd.read_table(filepath, sep="\t")
    df_sam["Algo"] = algo
    return df_sam.loc[:,["Sample_Id", "Algo"]]


def extract_commonly_analyzed_samples(df_sam, algos):
    df_sam_agg = df_sam.sort_values(by="Algo").groupby("Sample_Id").agg({"Algo": "|".join}).reset_index()
    sam_all = df_sam_agg.loc[df_sam_agg["Algo"]=="|".join(algos), "Sample_Id"].values.tolist()
    print("-INFO: there are %d sample analyzed by all of the following algorithms:" % len(sam_all))
    print("\t" + "\n\t".join(algos))
    return sam_all


def convert_type_to_str(x):
    try:
        y = "%d" % int(x)
    except:
        try:
            y = "%d" % float(x)
            if y=="nan":
                y = np.nan
        except:
            y = x
    return y


def add_id_with_breakpoint(df, col_id_w_brk):
    df_c = df.copy()
    df_c["Breakpoint_1"] = df_c["Breakpoint_1"].apply(convert_type_to_str).fillna("N/A")
    df_c["Breakpoint_2"] = df_c["Breakpoint_2"].apply(convert_type_to_str).fillna("N/A")
    df[col_id_w_brk] = df_c["Fusion_Id"] + "_" + df_c[["Breakpoint_1", "Breakpoint_2"]].apply("--".join, axis=1)
    return df


def add_aggregated_callers(df, col_id, cols_call, algos_to_prefix):
    dfs_algo = []
    for algo, prefix in algos_to_prefix.items():
        col_call_algo = prefix + "_" + "Call"
        if col_call_algo in cols_call:
            df_algo = df.loc[df[col_call_algo]==1][[col_id]]
            df_algo["Algo"] = algo
            dfs_algo.append(df_algo)
    df_algo = pd.concat(dfs_algo,axis=0)
    df_algo = df_algo.sort_values(by="Algo").groupby(col_id).agg({"Algo":"|".join}).reset_index()
    return df.merge(df_algo, how="left", on=col_id)


def select_combination(df, algos, n=None, n_exact=False):
    if n is not None:
        if n_exact:
            mask = df[algos].sum(axis=1)==n
        else:
            mask = df[algos].sum(axis=1)>=n
    else:
        mask = pd.Series(True, index=df.index)
        for algo in algos:
            mask = mask & (df[algo]==1)
    return mask


def main(args):
    algos_to_prefix = {"arriba": "AR",
                       "ericscript": "ES",
                       "fusioncatcher": "FC",
                       "pizzly": "PZ",
                       "squid": "SQ",
                       "starfusion": "SF",
                       "deepest_pnas_2019": "Deepest_Pnas_2019",
                       "prada_nar_2018": "Prada_Nar_2018",
                       "starfusion_cell_2018": "Starfusion_Cell_2018"}

    # load fusions
    df_fus = load_table(args.input, low_memory=False)

    # reset FFPM >= 0.1 for starfusion
    if "SF_Call" in df_fus.columns:
        df_fus.loc[(df_fus["SF_Call"]==1) & (df_fus["SF_FFPM"] < 0.1), "SF_Call"] = 0

    # choose algos
    if args.cohort=="tcga":
        algos = ["deepest_pnas_2019", "prada_nar_2018", "starfusion_cell_2018"]
    else:
        algos = ["arriba", "ericscript", "pizzly", "starfusion"]

    # load samples
    df_sam = pd.concat([load_table_samples(algo=algo, cohort=args.cohort) for algo in algos])

    # select data for commonly analyzed samples
    sam_all = extract_commonly_analyzed_samples(df_sam=df_sam, algos=algos)
    df_fus = df_fus.loc[df_fus["Sample_Id"].isin(sam_all)]

    # filter out potential false-positives and fusions involving
    # - uncharacterized genes
    # - immunoglobulin genes
    # - mitochondrial genes
    # - etc
    df_fus = filter_out_fusions(df_fus, tag_only=True)

    # filter on whitelist
    df_fus = filter_in_whitelist(df_fus, tag_only=True)

    # add ids with breakpoint
    col_id_wo_brk = "Fusion_Id"
    col_id_w_brk = "Fusion_Id_With_Breakpoint"
    df_fus = add_id_with_breakpoint(df=df_fus, col_id_w_brk=col_id_w_brk)

    # add row identifier and merge identifier
    col_row_id = "Row_Id"
    col_mer_id = "Merge_Id"
    df_fus[col_row_id] = df_fus[["Sample_Id", col_id_w_brk]].apply("_".join, axis=1)
    df_fus[col_mer_id] = df_fus[["Sample_Id", col_id_wo_brk]].apply("_".join, axis=1)

    # add aggregate callers per row
    if args.cohort=="tcga":
        cols_call = ["Deepest_Pnas_2019_Call", "Prada_Nar_2018_Call", "Starfusion_Cell_2018_Call"]
    else:
        cols_call = ["AR_Call", "ES_Call", "PZ_Call", "SF_Call"]
    df_fus = add_aggregated_callers(df_fus, col_id=col_row_id, cols_call=cols_call, algos_to_prefix=algos_to_prefix)

    # build table of fusions without breakpoints for filtering
    # if col_id does not contain breakpoint information, it may happen that for one gene pair
    # in one sample, a breakpoint was called by a given caller but not other breakpoints. In order to
    # harmonize calls at the level of col_id, if any breakpoint is identified for this gene pair
    # by the caller, the gene pair will be considered as called by the caller.
    df_fus_ids = df_fus[[col_mer_id]+cols_call].drop_duplicates()
    df_fus_ids = df_fus_ids.drop_duplicates().groupby(col_mer_id).sum()
    df_fus_ids = (df_fus_ids > 0).astype(int).reset_index().copy()
    df_fus_ids = df_fus_ids.loc[df_fus_ids.loc[:,cols_call].sum(axis=1)>0].copy()
    df_fus_ids = add_aggregated_callers(df_fus_ids, col_id=col_mer_id, cols_call=cols_call, algos_to_prefix=algos_to_prefix)

    assert df_fus_ids.shape[0]==df_fus_ids[col_mer_id].nunique()

    # filter on algorithms
    if args.cohort=="tcga":
        mask_1 = select_combination(df_fus_ids, algos=["Starfusion_Cell_2018_Call"])
        mask_2 = select_combination(df_fus_ids, algos=["Deepest_Pnas_2019_Call", "Prada_Nar_2018_Call"])
        mask = mask_1 | mask_2
    else:
        mask_1 = select_combination(df_fus_ids, algos=["AR_Call", "ES_Call"])
        mask_2 = select_combination(df_fus_ids, algos=["PZ_Call", "SF_Call"])
        mask = mask_1 | mask_2

    tag = "Algorithms_Combination"
    df_fus_ids.loc[~mask, "FILTER_IDS"] = tag

    # filter on table with breakpoints
    df_fus = df_fus.rename(columns={"Algo": "Algo_Wt_Breakpoint"})
    df_fus_ids = df_fus_ids.rename(columns={"Algo": "Algo_Wo_Breakpoint"})
    df_fus = df_fus.merge(df_fus_ids[[col_mer_id, "Algo_Wo_Breakpoint", "FILTER_IDS"]], how="left", on=col_mer_id)

    # merge FILTER, FILTER_IDS
    mask_null = df_fus["FILTER"].isnull()
    mask_ids = ~df_fus["FILTER_IDS"].isnull()
    tag = "Algorithms_Combination"
    df_fus.loc[(mask_ids) & (~mask_null), "FILTER"] += ";%s" % tag
    df_fus.loc[(mask_ids) & (mask_null), "FILTER"] = tag
    del df_fus["FILTER_IDS"]
    del df_fus["Merge_Id"]

    # filter on samples in design
    if args.cohort=="tcga_validation":
        filepath_bio = "../../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv"
    else:
        filepath_bio = "../../../data/%s/clinical/curated/bio_%s_in_design_curated.tsv" % (args.cohort, args.cohort)
    df_bio = pd.read_table(filepath_bio)
    df_bio = df_bio.loc[df_bio["Sample_Type"].isin(["RNA_N", "RNA_T"])]

    mask = df_fus["Sample_Id"].isin(df_bio["Sample_Id"])
    df_fus = df_fus.loc[mask]
    print("-INFO: %d/%d lines are selected from samples present in %s" % (sum(mask), mask.shape[0], filepath_bio))

    # remove per-caller columns
    df_fus = df_fus[[x for x in df_fus if not x.endswith("_Call")]]

    # move id columns to beginning
    cols_ids = [x for x in df_fus if x.endswith("_Id") or "Fusion_Id" in x]
    cols_oth = [x for x in df_fus if x not in cols_ids]
    df_fus = df_fus[cols_ids+cols_oth]
    df_fus = df_fus.sort_values(by="Row_Id")

    # make table of fusions filters
    cols_filters = ["Sample_Id", "Fusion_Id", "Fusion_Id_With_Breakpoint", "FILTER"]
    df_fus.loc[df_fus["FILTER"].isnull(), "FILTER"] = "PASS"
    df_fus_filters = df_fus[cols_filters].copy()
    df_fus_filters.to_csv(args.output_fus_filters, sep="\t", index=False)
    print("-table of fusions filters saved at %s" % args.output_fus)

    # select only PASS
    df_fus = df_fus.loc[df_fus["FILTER"]=="PASS"].copy()

    # format some numeric columns
    cols = ["Breakpoint_1", "Breakpoint_2"] + [x for x in df_fus if x.endswith("_Reads")]
    for col in cols:
        df_fus[col] = df_fus[col].apply(convert_type_to_str)

    # save fusions
    df_fus.to_csv(args.output_fus, sep="\t", index=False)
    print("-table of filtered fusions saved at %s" % args.output_fus)

    # save table sample_list
    sam_all_bio = [x for x in sam_all if x in df_bio["Sample_Id"].unique()]
    pd.DataFrame({"Sample_Id": sam_all_bio}).to_csv(args.output_sam, sep="\t", index=False)
    print("-table of samples saved at %s" % args.output_sam)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a table of filtered fusion events.')
    parser.add_argument('--cohort', type=str, help='Name of the cohort.', default="met500")
    parser.add_argument('--input', type=str, help='Path to annotated fusion calls.',
                        default="../../../data/met500/rna/fusions/met500_annotated.tsv.gz")
    parser.add_argument('--output_sam', type=str, help='Path to output table of samples.',
           default="../../../data/met500/rna/fusions/sample_list.tsv")
    parser.add_argument('--output_fus_filters', type=str, help='Path to output table of scores.',
           default="../../../data/met500/rna/fusions/met500_filters.tsv.gz")
    parser.add_argument('--output_fus', type=str, help='Path to output table of scores.',
           default="../../../data/met500/rna/fusions/met500_annotated_filtered.tsv.gz")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")

    main(args)
