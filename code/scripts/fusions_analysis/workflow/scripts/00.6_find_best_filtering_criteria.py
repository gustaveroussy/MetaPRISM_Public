# -*- coding: utf-8 -*-
"""
@created: Sep 10 2021
@modified: Feb 12 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France
"""

import argparse
import itertools
from joblib import Parallel, delayed
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
from utils import load_table_samples, extract_commonly_analyzed_samples, add_aggregated_callers
from pyprism.data import load_bio, load_table

# functions ============================================================================================================

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


def add_high_confidence_calls_per_caller(df_a, df_b, cols_call_a, col_id_w_brk):
    col_row_id_hc = "Row_Id_HC"
    df_a[col_row_id_hc] = df_a[["Sample_Id", col_id_w_brk]].apply("_".join, axis=1)
    df_b[col_row_id_hc] = df_b[["Sample_Id", col_id_w_brk]].apply("_".join, axis=1)

    # high-confidence calls are defined as calls in true set that are called by all callers
    mask_a_hc = df_a[cols_call_a].sum(axis=1)==len(cols_call_a)
    ids_a_hc = df_a.loc[mask_a_hc, col_row_id_hc]

    cols_call_hc_b = []

    # arriba 
    ar_hc_mask = df_b["AR_Confidence"]=="high"
    df_b["AR_HC_Call"] = 0
    df_b.loc[ar_hc_mask, "AR_HC_Call"] = 1
    cols_call_hc_b.append("AR_HC_Call")

    # ericscript
    es_score_hc_threshold = df_b.loc[df_b[col_row_id_hc].isin(ids_a_hc)]["ES_Score"].mean()
    es_hc_mask = df_b["ES_Score"]>=es_score_hc_threshold
    df_b["ES_HC_Call"] = 0
    df_b.loc[es_hc_mask, "ES_HC_Call"] = 1
    cols_call_hc_b.append("ES_HC_Call")

    # starfusion
    sf_ffpm_hc_threshold = df_b.loc[df_b[col_row_id_hc].isin(ids_a_hc)]["SF_FFPM"].median()
    sf_hc_mask = df_b["SF_FFPM"]>=sf_ffpm_hc_threshold
    df_b["SF_HC_Call"] = 0
    df_b.loc[sf_hc_mask, "SF_HC_Call"] = 1
    cols_call_hc_b.append("SF_HC_Call")

    del df_a[col_row_id_hc]
    del df_b[col_row_id_hc]

    return df_b, cols_call_hc_b



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


def get_all_selections(df, cols_call, suffix=""):
    selections = {}
    for i in range(1, len(cols_call)+1):
        # all combination of callers. For n callers, there are 2^n -1 possible combinations
        for subset in itertools.combinations(cols_call, i):
            selections[" & ".join(sorted(list(subset)))] = select_combination(df, subset)


        # calls by i callers or more. For n callers, there are n possible combinations
        if i < len(cols_call):
            selections["%d or more%s" % (i, suffix)] = select_combination(df, cols_call, n=i, n_exact=False)

        # all combinations of callers except all callers and/or all combinations of callers not in first combination
        if i < len(cols_call):
            for subset in itertools.combinations(cols_call, i):
                mask_first = select_combination(df, subset)

                cols_call_left = list(set(cols_call).difference(set(subset)))
                for i_left in range(1, len(cols_call_left)+1):
                    for subset_left in itertools.combinations(cols_call_left, i_left):
                        mask_left = select_combination(df, subset_left)

                        name = " & ".join(sorted(list(subset))) + " || " + " & ".join(sorted(list(subset_left)))
                        mask = mask_first | mask_left

                        selections[name] = mask

                    if i_left > 1 and i_left < len(cols_call_left):
                        # calls by i_left callers or more
                        mask_left = select_combination(df, cols_call_left, n=i_left, n_exact=False)

                        # union of calls by subset callers and calls by i_left callers or more
                        name = " & ".join(sorted(list(subset))) + " || " + "%d or more%s" % (i_left, suffix)
                        mask = mask_first | mask_left
                        selections[name] = mask

                        # intersection of calls by subset callers and calls by i_left callers or more
                        name = " & ".join(sorted(list(subset))) + " && " + "%d or more%s" % (i_left, suffix)
                        mask = mask_first & mask_left
                        selections[name] = mask
    return selections


def get_scores_combination(df_a,df_b,col_id,a_is_truth=False):
    A = df_a[col_id]
    B = df_b[col_id]
    if a_is_truth:
        X,Y=set(A),set(B)
    else:
        X,Y=set(B),set(A)

    scores = {}
    scores["N_True"] = len(X)
    scores["N_Pred"] = len(Y)
    scores["Prop_True_Found"] = len(X.intersection(Y))/len(X)
    scores["Prop_False_Found"] = len(Y.difference(X))/len(X)
    scores["DSC"] = 2*len(set(X).intersection(set(Y)))/(len(set(X)) + len(set(Y)))

    return scores



def get_scores_combination_fixed_ab(df_a_ids, df_b_ids, col_row_id, name_a, name_b, mask_a, mask_b, selections_b_hc):
    scores_all = []
    scores = get_scores_combination(df_a_ids.loc[mask_a],df_b_ids.loc[mask_b],col_id=col_row_id,a_is_truth=True)
    scores["Name_A"] = name_a
    scores["Name_B"] = name_b
    scores_all.append(scores)

    for name_b_hc, mask_b_hc in selections_b_hc.items():
        scores = get_scores_combination(df_a_ids.loc[mask_a],df_b_ids.loc[mask_b | mask_b_hc],
                                        col_id=col_row_id,a_is_truth=True)
        scores["Name_A"] = name_a
        scores["Name_B"] = name_b + " || " + name_b_hc
        scores_all.append(scores)

    return scores_all



def detailed_comparison(df_a_ids, df_b_ids, df_a, df_b, col_row_id, name_a, name_b, mask_a, mask_b):
    cols_b_not_a = [x for x in df_b if x not in df_a]
    cols_a_not_b = [x for x in df_a if x not in df_b]

    df_a_ids_for_merge = df_a_ids[[col_row_id, "Algo"]].rename(columns={"Algo": "Algo_A"})
    df_a_ids_for_merge = df_a_ids_for_merge.merge(df_a, how="left", on=col_row_id)
    df_b_ids_for_merge = df_b_ids[[col_row_id, "Algo"]].rename(columns={"Algo": "Algo_B"})
    df_b_ids_for_merge = df_b_ids_for_merge.merge(df_b, how="left", on=col_row_id)

    df_a_ids_sel = df_a_ids.loc[mask_a]
    df_b_ids_sel = df_b_ids.loc[mask_b]
    ids_a_and_b = set(df_a_ids_sel[col_row_id]).intersection(set(df_b_ids_sel[col_row_id]))
    ids_a_not_b = set(df_a_ids_sel[col_row_id]).difference(set(df_b_ids_sel[col_row_id]))
    ids_b_not_a = set(df_b_ids_sel[col_row_id]).difference(set(df_a_ids_sel[col_row_id]))

    df_a_ids_a_and_b = df_a_ids_sel.loc[df_a_ids_sel[col_row_id].isin(list(ids_a_and_b))][[col_row_id]]
    df_a_ids_a_and_b = df_a_ids_a_and_b.merge(df_a_ids_for_merge, how="left", on=col_row_id)
    df_a_ids_a_and_b = df_a_ids_a_and_b.merge(df_b_ids_for_merge[[col_row_id]+cols_b_not_a], how="left", on=col_row_id)

    df_a_ids_a_not_b = df_a_ids_sel.loc[df_a_ids_sel[col_row_id].isin(list(ids_a_not_b))][[col_row_id]]
    df_a_ids_a_not_b = df_a_ids_a_not_b.merge(df_a_ids_for_merge, how="left", on=col_row_id)
    df_a_ids_a_not_b = df_a_ids_a_not_b.merge(df_b_ids_for_merge[[col_row_id]+cols_b_not_a], how="left", on=col_row_id)

    df_b_ids_b_not_a = df_b_ids_sel.loc[df_b_ids_sel[col_row_id].isin(list(ids_b_not_a))][[col_row_id]]
    df_b_ids_b_not_a = df_b_ids_b_not_a.merge(df_b_ids_for_merge, how="left", on=col_row_id)
    df_b_ids_b_not_a = df_b_ids_b_not_a.merge(df_a_ids_for_merge[[col_row_id]+cols_a_not_b], how="left", on=col_row_id)

    # save tables with calls in a not b and vice-versa
    folder = "../../../results/fusions_analysis/validation/%s_vs_%s" % (name_a, name_b)
    os.makedirs(folder, exist_ok=True)
    df_a_ids_a_not_b.to_excel(os.path.join(folder, "a_ids_a_not_b.xlsx"), index=False)
    df_b_ids_b_not_a.to_excel(os.path.join(folder, "b_ids_b_not_a.xlsx"), index=False)

    # draw plots comparing metrics distribution between true positives, false-positives and false-negatives
    # false positives
    df_tp = df_a_ids_a_and_b
    df_fp = df_b_ids_b_not_a
    df_fn = df_a_ids_a_not_b

    for col in ["ES_Score", "SF_FFPM"]:
        colors = {"TP": "#7BE0AD", "FP": "#FAAA8D", "FN": "#7F7CAF"}

        y_tp = df_tp[col].dropna().to_frame(col)
        y_tp["Label"] = "TP"
        y_fn = df_fn[col].dropna().to_frame(col)
        y_fn["Label"] = "FN"
        y_fp = df_fp[col].dropna().to_frame(col)
        y_fp["Label"] = "FP"
        y_all = pd.concat((y_tp, y_fn, y_fp), axis=0).reset_index(drop=True)
        y_all[col] = np.log10(y_all[col])

        fig, ax = plt.subplots(figsize=(12,8))
        sns.kdeplot(data=y_all, x=col, hue="Label", cut=0, fill=True, common_norm=False, alpha=0.7)
        fig.canvas.draw()
        locs, labels = plt.xticks()
        # Or for scientific notation:
        ax.set(xticklabels=["$10^{" + i.get_text() + "}$" for i in labels])
        plt.savefig(os.path.join(folder, "kdeplot_%s_tp_fp_fn.pdf" % col))

    for prefix in ["AR", "ES", "SF"]:
        cols = ["%s_Spanning_Reads" % prefix, "%s_Split_Reads" % prefix]

        y_tp = df_tp[cols].dropna()
        y_tp["Label"] = "TP"
        y_fn = df_fn[cols].dropna()
        y_fn["Label"] = "FN"
        y_fp = df_fp[cols].dropna()
        y_fp["Label"] = "FP"
        y_all = pd.concat((y_tp, y_fn, y_fp), axis=0).reset_index(drop=True)
        for col in cols:
            y_all[col] = np.log10(y_all[col]+1)


        fig,ax = plt.subplots(1,1,figsize=(8,8))
        for label in y_all["Label"].unique():
            mask = y_all["Label"]==label
            ax.scatter(y_all.loc[mask, cols[0]], y_all.loc[mask, cols[1]], color=colors[label], s=25, label=label)
        ax.set_xlabel(cols[0], fontsize=16)
        ax.set_ylabel(cols[0], fontsize=16)
        ax.legend(loc="best", fontsize=16)
        fig.canvas.draw()
        locs, labels = plt.xticks()
        # Or for scientific notation:
        ax.set(xticklabels=["$10^{" + i.get_text() + "}$" for i in labels])
        locs, labels = plt.yticks()
        # Or for scientific notation:
        ax.set(yticklabels=["$10^{" + i.get_text() + "}$" for i in labels])
        plt.savefig(os.path.join(folder, "scatter_%s_%s_tp_fp_fn.pdf" % (cols[0], cols[1])))


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

    # load
    cohort_a = "tcga"
    cohort_b = "tcga_validation"
    algos_a = args.algos_tcga
    algos_b = args.algos_tcga_validation

    # load fusions
    df_a = load_table(args.fusions_tcga, low_memory=False)
    df_b = load_table(args.fusions_tcga_validation, low_memory=False)

    # reset FFPM >= 0.1 for starfusion
    df_b.loc[(df_b["SF_Call"]==1) & (df_b["SF_FFPM"] < 0.1), "SF_Call"] = 0

    # load samples
    df_sam_a = pd.concat([load_table_samples(algo=algo, cohort=cohort_a) for algo in algos_a])
    df_sam_b = pd.concat([load_table_samples(algo=algo, cohort=cohort_b) for algo in algos_b])
    print("-INFO: %d unique samples in %s" % (df_sam_a["Sample_Id"].nunique(), cohort_a))
    print("-INFO: %d unique samples in %s" % (df_sam_b["Sample_Id"].nunique(), cohort_b))

    # select data for commonly analyzed samples
    sam_com = extract_commonly_analyzed_samples(df_sam_a=df_sam_a, df_sam_b=df_sam_b)
    print("-INFO: %d commonly analyzed samples btw %s and %s" % (len(sam_com), cohort_a, cohort_b))
    df_a = df_a.loc[df_a["Sample_Id"].isin(sam_com)].copy()
    df_b = df_b.loc[df_b["Sample_Id"].isin(sam_com)].copy()

    # filter out potential false-positives and fusions involving
    # - uncharacterized genes
    # - immunoglobulin genes
    # - mitochondrial genes
    # - etc
    df_a = filter_out_fusions(df_a)
    df_b = filter_out_fusions(df_b)

    # filter on whitelist
    if args.use_whitelist.lower()=="yes":
        df_a = filter_in_whitelist(df_a)
        df_b = filter_in_whitelist(df_b)

    # choose how a fusion is identified
    # add ids with breakpoint
    col_id_wo_brk = "Fusion_Id"
    col_id_w_brk = "Fusion_Id_With_Breakpoint"
    df_a = add_id_with_breakpoint(df=df_a, col_id_w_brk=col_id_w_brk)
    df_b = add_id_with_breakpoint(df=df_b, col_id_w_brk=col_id_w_brk)

    if args.use_breakpoints.lower()=="yes":
        col_id = col_id_w_brk
    else:
        col_id = col_id_wo_brk

    # add row identifier
    col_row_id = "Row_Id"
    df_a[col_row_id] = df_a[["Sample_Id", col_id]].apply("_".join, axis=1)
    df_b[col_row_id] = df_b[["Sample_Id", col_id]].apply("_".join, axis=1)

    cols_call_a = ["%s_Call" % algos_to_prefix[algo] for algo in algos_a]
    cols_call_b = ["%s_Call" % algos_to_prefix[algo] for algo in algos_b]

    # determine high-confidence thresholds for ES_Score, SF_FFPM so that fusions called only by one of this
    # algorithm at a high-confidence will be rescued
    df_b, cols_call_hc_b = add_high_confidence_calls_per_caller(df_a, df_b, cols_call_a, col_id_w_brk)

    # if col_id does not contain breakpoint information, it may happen that for one gene pair
    # in one sample, a breakpoint was called by a given caller but not other breakpoints. In order to
    # harmonize calls at the level of col_id, if any breakpoint is identified for this gene pair
    # by the caller, the gene pair will be considered as called by the caller.
    df_a_ids = df_a[[col_row_id]+cols_call_a].drop_duplicates()
    df_a_ids = df_a_ids.drop_duplicates().groupby(col_row_id).sum()
    df_a_ids = (df_a_ids > 0).astype(int).reset_index().copy()
    df_a_ids = df_a_ids.loc[df_a_ids.loc[:,cols_call_a].sum(axis=1)>0].copy()
    df_b_ids = df_b[[col_row_id]+cols_call_b+cols_call_hc_b].drop_duplicates()
    df_b_ids = df_b_ids.drop_duplicates().groupby(col_row_id).sum()
    df_b_ids = (df_b_ids > 0).astype(int).reset_index().copy()
    df_b_ids = df_b_ids.loc[df_b_ids.loc[:,cols_call_b+cols_call_hc_b].sum(axis=1)>0].copy()

    assert df_a_ids.shape[0]==df_a_ids[col_row_id].nunique()
    assert df_b_ids.shape[0]==df_b_ids[col_row_id].nunique()

    # add aggregate callers per row
    df_a_ids = add_aggregated_callers(df_a_ids, col_id=col_row_id, cols_call=cols_call_a,
                                      algos_to_prefix=algos_to_prefix)
    df_b_ids = add_aggregated_callers(df_b_ids, col_id=col_row_id, cols_call=cols_call_b,
                                      algos_to_prefix=algos_to_prefix)

    # info about ids in common and different
    ids_a_and_b = set(df_a_ids[col_row_id]).intersection(set(df_b_ids[col_row_id]))
    ids_a_not_b = set(df_a_ids[col_row_id]).difference(set(df_b_ids[col_row_id]))
    ids_b_not_a = set(df_b_ids[col_row_id]).difference(set(df_a_ids[col_row_id]))
    print("-INFO: %d row in %s and %s" % (len(ids_a_and_b), cohort_a, cohort_b))

    print("-INFO: %d row in %s not in %s" % (len(ids_a_not_b), cohort_a, cohort_b))
    cnt_a_not_b = df_a_ids.loc[df_a_ids[col_row_id].isin(ids_a_not_b), "Algo"].value_counts()
    print("--for ids in %s not in %s" % (cohort_a, cohort_b))
    for algo, cnt in zip(cnt_a_not_b.index, cnt_a_not_b.values):
        print("\t%d ids called by %s" % (cnt, algo))

    print("-INFO: %d row in %s not in %s" % (len(ids_b_not_a), cohort_b, cohort_a))
    cnt_b_not_a = df_b_ids.loc[df_b_ids[col_row_id].isin(ids_b_not_a), "Algo"].value_counts()
    print("--for ids in %s not in %s" % (cohort_b, cohort_a))
    for algo, cnt in zip(cnt_b_not_a.index, cnt_b_not_a.values):
        print("\t%d ids called by %s" % (cnt, algo))

    # find best combinations of callers/filters
    selections_a = get_all_selections(df=df_a_ids, cols_call=cols_call_a)
    selections_b = get_all_selections(df=df_b_ids, cols_call=cols_call_b)
    selections_b_hc = get_all_selections(df=df_b_ids, cols_call=cols_call_hc_b, suffix=" HC")
    df_scores = pd.DataFrame(columns=["Name_A", "Name_B", "DSC", "Prop_True_Found", "Prop_False_Found", "N_True",
                                      "N_Pred"])

    for name_a, mask_a in tqdm(selections_a.items()):
        task_parallel = lambda name_b: get_scores_combination_fixed_ab(df_a_ids, df_b_ids, col_row_id, name_a,
                                                                       name_b, mask_a, selections_b[name_b],
                                                                       selections_b_hc)
        scores_all = Parallel(n_jobs=args.n_jobs)(delayed(task_parallel)(name_b) for name_b in selections_b.keys())

        for scores in scores_all:
            df_scores = df_scores.append(scores, ignore_index=True)

    # save scores
    df_scores = df_scores.sort_values(by="DSC", ascending=False)
    df_scores.to_csv(args.output_table, sep="\t", index=False, float_format="%.2f")
    print("-table of scores saved at %s" % args.output_table)

    # # analyze in details overlap between 2 selections
    # name_a = 'Deepest_Pnas_2019_Call & Prada_Nar_2018_Call || Starfusion_Cell_2018_Call'
    # mask_a = selections_a[name_a]
    # name_b = 'AR_Call & ES_Call || PZ_Call & SF_Call'
    # mask_b = selections_b[name_b]
    # detailed_comparison(df_a_ids, df_b_ids, df_a, df_b, col_row_id, name_a, name_b, mask_a, mask_b)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make a curated table of the samples in the study design.')
    parser.add_argument('--fusions_tcga', type=str, help='Path to annotated fusion calls for tcga.',
                        default="../../../data/tcga/rna/fusions/tcga_annotated.tsv.gz")
    parser.add_argument('--fusions_tcga_validation', type=str, help='Path to annotated fusion calls for tcga validation.',
                        default="../../../data/tcga_validation/rna/fusions/tcga_validation_annotated.tsv.gz")
    parser.add_argument("--algos_tcga", type=str, nargs="+", help="Algos analyzed for tcga data",
                      default=["deepest_pnas_2019", "prada_nar_2018", "starfusion_cell_2018"])
    parser.add_argument("--algos_tcga_validation", type=str, nargs="+", help="Algos analyzed for tcga validation data",
                      default=["arriba", "ericscript", "pizzly", "starfusion"])
    parser.add_argument('--use_whitelist', type=str, default="Yes",
        help='If set to "Yes", only fusions listed in the the whitelist table will be considered')
    parser.add_argument('--use_breakpoints', type=str, default="Yes",
        help='If set to "Yes", fusion call is identical if partners and breakpoints are identical.')
    parser.add_argument('--n_jobs', type=int, default=4, help='Number of cores available for parallel computations')
    parser.add_argument('--output_table', type=str, help='Path to output table of scores .',
           default="../../../results/fusions_analysis/validation/table_scores_use_whitelist_yes_use_breakpoints_no.tsv")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")

    main(args)
