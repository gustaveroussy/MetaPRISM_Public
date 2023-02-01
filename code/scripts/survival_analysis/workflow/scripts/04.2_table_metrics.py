# -*- coding: utf-8 -*-
"""
@created: May 17 2022
@modified: Aug 08 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Save a summary table of c-indexes.
"""

import argparse
import numpy as np
import pandas as pd
import os
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import sys
sys.path.append("workflow/functions")

from utils import load_met_cov, load_cov_ori

def select_rows_met(df_met, time_horizon):
    df_met_rep = df_met.copy()

    # choose time horizon for metrics
    if time_horizon not in df_met_rep["Time"].unique().tolist():
        preselected_times = df_met_rep["Time"].unique().tolist()
        raise ValueError("-the specified time horizon %d is not in the preselected list of time horizons:\n\t%s" %
                         (time_horizon, preselected_times))
    mask_time = df_met_rep["Time"]==time_horizon
    print("INFO: selecting %d/%d lines for time horizon %d." % (sum(~mask_time), len(mask_time), time_horizon))

    df_met_rep = df_met_rep.loc[mask_time].copy()

    # choose only test splits
    mask_test = df_met_rep["Split"]=="Test"
    print("INFO: selecting %d/%d lines for split Test." % (sum(~mask_test), len(mask_test)))
    df_met_rep = df_met_rep.loc[mask_test].copy()

    return df_met_rep


def main(args):
    # choose default values for models and selection unless values were specified by user
    if len(args.names_samples)==0:
        samples = ["BLCA_dna_and_rna", "BRCA_dna_and_rna", "LUAD_dna_and_rna", "PAAD_dna_and_rna", "PRAD_dna_and_rna",
                   "dna_cohort_dna_and_rna"]
    else:
        samples = args.names_samples

    if len(args.names_features)==0:
        features = ["cln_grim_0", "cln_biol_0", "cln_biol_1", "cln_biol_2",
                    "cln_biol_1_dna_1", "cln_biol_1_dna_2", "cln_biol_1_dna_3", "cln_biol_1_dna_4",
                    "cln_biol_1_dna_5", "cln_biol_1_rna_1", "cln_biol_1_rna_2",
                    "cln_biol_1_dna_1_rna_1", "cln_biol_1_dna_5_rna_1", "cln_dna_1", "cln_rna_1", "cln_dna_1_rna_1"]
    else:
        features = args.names_features

    if len(args.names_models)==0:
        models = ["coxph_standard", "coxph_ridge", "coxph_lasso"]
    else:
        models = args.names_models

    if len(args.names_selections)==0:
        selections = ["none", "lasso", "10pct_fdr"]
    else:
        selections = args.names_selections

    # load pooled metrics and coefficients estimates
    dfs_met = []
    dfs_pval = []
    dfs_cov_fin = []
    dfs_cov_ori = []

    for sample in samples:
        dir_data = os.path.join(args.dir_results, "data_%s/sub_features" % args.name_cohort, sample)
        dir_models = os.path.join(args.dir_results, "models_%s" % args.name_cohort, sample)

        df_met, df_cov_fin = load_met_cov(dir_models=dir_models,
                                          dir_data=dir_data,
                                          names_features=features,
                                          names_models=models,
                                          names_selections=selections,
                                          pool_level="rep")
        df_met = select_rows_met(df_met=df_met, time_horizon=args.time_horizon)


        # load across repeats in order to compute pvalues from wilcoxon tests
        df_met_all, _ = load_met_cov(dir_models=dir_models,
                                     dir_data=dir_data,
                                     names_features=features,
                                     names_models=models,
                                     names_selections=selections,
                                     pool_level="imp")
        df_met_all = select_rows_met(df_met=df_met_all, time_horizon=args.time_horizon)

        mask_f = df_met_all["Features"]=="cln_biol_1"
        mask_s = df_met_all["Selection"]=="none"
        mask_m = df_met_all["Model"]=="coxph_standard"
        mask_a = mask_f & mask_s & mask_m
        mets_a = df_met_all.loc[mask_a]["C_ipcw"].values

        df_pval = pd.DataFrame(columns=["Features", "Selection", "Model", "pval"])
        df_comb = df_met_all[["Features", "Selection", "Model"]].drop_duplicates()
        for _, comb in df_comb.iterrows():
            if comb["Features"]=="cln_grim_0" or comb["Features"]=="cln_biol_1":
                pass
            else:
                mask_f = df_met_all["Features"]==comb["Features"]
                mask_s = df_met_all["Selection"]==comb["Selection"]
                mask_m = df_met_all["Model"]==comb["Model"]
                mask_b = mask_f * mask_s * mask_m

                mets_b = df_met_all.loc[mask_b]["C_ipcw"].values
                comb["pval"] = wilcoxon(x=mets_b, y=mets_a, alternative="greater")[-1]
                df_pval = df_pval.append(comb, ignore_index=True)


        df_cov_ori = load_cov_ori(dir_data=dir_data,
                              names_features=features)

        df_met["Sample"] = sample
        df_pval["Sample"] = sample
        df_cov_fin["Sample"] = sample
        df_cov_ori["Sample"] = sample

        dfs_met.append(df_met)
        dfs_pval.append(df_pval)
        dfs_cov_fin.append(df_cov_fin)
        dfs_cov_ori.append(df_cov_ori)

    # aggregate
    df_met = pd.concat(dfs_met)
    df_pval = pd.concat(dfs_pval)
    df_cov_ori = pd.concat(dfs_cov_ori)
    df_cov_fin = pd.concat(dfs_cov_fin)

    # remove Survival_Time and Survival_Status
    covs_rmv = ["Survival_Time", "Survival_Status"]
    df_cov_ori = df_cov_ori.loc[~df_cov_ori["Covariate"].isin(covs_rmv)].copy()
    df_cov_fin = df_cov_fin.loc[~df_cov_fin["Covariate"].isin(covs_rmv)].copy()

    # column group by
    cols_gby_met = ["Sample", "Features", "Selection", "Model"]
    cols_gby_cov_1 = ["Sample", "Features"]
    cols_gby_cov_2 = ["Sample", "Features", "Selection"]

    df_met_agg = df_met[cols_gby_met + ["C_ipcw"]].copy()

    # add Nb Features A, Nb Features B, Nb Features Selected 20% and Data Type
    name_cov_ori = "Nb Features A"
    name_cov_1 = "Nb Features B"
    name_cov_2 = "Nb Features Selected %d%%" % (args.threshold*100)

    df_cov_ori_agg = df_cov_ori.groupby(cols_gby_cov_1).size().to_frame(name_cov_ori).reset_index()
    df_cov_fin_1 = df_cov_fin[["Sample", "Features", "Covariate"]].drop_duplicates()
    df_cov_fin["Selected_WE"] = df_cov_fin[["N_Selected", "N_Warning", "N_Error"]].sum(axis=1)/df_cov_fin["N_Repeats"]
    df_cov_fin_sel = df_cov_fin.loc[df_cov_fin["Selected_WE"] >= args.threshold]
    df_cov_fin_2 = df_cov_fin_sel[["Sample", "Features", "Selection", "Covariate"]].drop_duplicates()

    df_cov_fin_1_agg = df_cov_fin_1.groupby(cols_gby_cov_1).size().to_frame(name_cov_1).reset_index()
    df_cov_fin_2_agg = df_cov_fin_2.groupby(cols_gby_cov_2).size().to_frame(name_cov_2).reset_index()
    df_cov_fin_agg = df_cov_fin_1_agg.merge(df_cov_fin_2_agg, how="outer", on=cols_gby_cov_1)

    df_cov_agg  = df_cov_ori_agg.merge(df_cov_fin_agg, how="outer", on=cols_gby_cov_1)
    df_met_agg = df_met_agg.merge(df_cov_agg, how="left", on=cols_gby_cov_2)

    # add pval
    df_met_agg = df_met_agg.merge(df_pval, how="left", on=["Features", "Selection", "Model", "Sample"])
    df_met_agg["pval"] = multipletests(df_met_agg["pval"].fillna(1).values.tolist(), alpha=0.05, method="fdr_bh")[1]
    df_met_agg["C_ipcw"] = df_met_agg[["C_ipcw", "pval"]].apply(lambda x: "%.3f^*" % x["C_ipcw"] if x["pval"] < 0.05 else
                                                                "%.3f" % x["C_ipcw"], axis=1)

    # select columns and reformat
    df_met_agg = df_met_agg.rename(columns={"Nb Features A": "Nb Features"})
    cols_keep = ["Sample", "Features", "Nb Features", "Selection", "Model", "C_ipcw"]
    df_met_agg["Sample"] = df_met_agg["Sample"].apply(lambda x: x.split("_dna_and_rna")[0])

    # table full
    df_met_agg_full = df_met_agg[cols_keep].set_index(cols_keep[:-1]).unstack(level=0)
    df_met_agg_full.to_excel(args.output_full)
    print("-table saved at %s" % args.output_full)

    # table comb
    features_comb = ["cln_grim_0", "cln_biol_1", "cln_biol_1_dna_1", "cln_biol_1_dna_2", "cln_biol_1_rna_1",
                     "cln_biol_1_rna_2", "cln_biol_1_dna_1_rna_1", "cln_dna_1_rna_1"]
    selections_comb = ["lasso", "none"]
    models_comb = ["coxph_standard", "coxph_lasso"]

    mask_f = df_met_agg["Features"].isin(features_comb)
    mask_s = df_met_agg["Selection"].isin(selections_comb)
    mask_m = df_met_agg["Model"].isin(models_comb)
    mask = mask_f & mask_s & mask_m
    df_met_agg_comb = df_met_agg[mask].copy()

    # for grim 0 and biol 1, keep only "none" and "coxph_standard"
    features = ["cln_grim_0", "cln_biol_1", "cln_biol_1_dna_1", "cln_biol_1_rna_1", "cln_biol_1_dna_1_rna_1",
                "cln_dna_1_rna_1"]
    selections_keep = ["none"]
    models_keep = ["coxph_standard"]
    for feature in features:
        mask_f = df_met_agg_comb["Features"]==feature
        mask_s = df_met_agg_comb["Selection"].isin(selections_keep)
        mask_m = df_met_agg_comb["Model"].isin(models_keep)
        mask_bad = mask_f & (~mask_s | ~mask_m)
        df_met_agg_comb = df_met_agg_comb[~mask_bad]

    # for biol 1 dna 2, dont consider none
    features = ["cln_biol_1_dna_2", "cln_biol_1_rna_2"]
    selections_keep = ["lasso"]
    models_keep = ["coxph_standard"]
    for feature in features:
        mask_f = df_met_agg_comb["Features"]==feature
        mask_s = df_met_agg_comb["Selection"].isin(selections_keep)
        mask_m = df_met_agg_comb["Model"].isin(models_keep)
        mask_bad = mask_f & (~mask_s | ~mask_m)
        df_met_agg_comb = df_met_agg_comb[~mask_bad]

    df_met_agg_comb["Model"] = df_met_agg_comb["Model"].apply(lambda x: x.replace("coxph_", "Cox "))
    df_met_agg_comb["Features"] = df_met_agg_comb["Features"].replace({"cln_grim_0": "M1: GRIM disc.",
                                                                       "cln_biol_1": "M2: GRIM cont. + clinic",
                                                                       "cln_biol_1_dna_1": "M3: M2 + WES summ.",
                                                                       "cln_biol_1_dna_2": "M4: M3 + WES alt.",
                                                                       "cln_biol_1_rna_1": "M5: M2 + RNA summ.",
                                                                       "cln_biol_1_rna_2": "M6: M2 + RNA pthw.",
                                                                       "cln_biol_1_dna_1_rna_1": "M7: M3 + M5",
                                                                       "cln_dna_1_rna_1": "M7bis: M3 + M5 - GRIM cont"})

    df_met_agg_comb = df_met_agg_comb[cols_keep].set_index(cols_keep[:-1]).unstack(level=0)
    df_met_agg_comb.columns = df_met_agg_comb.columns.get_level_values(-1)
    df_met_agg_comb = df_met_agg_comb.rename(columns={"dna_cohort": "All"})
    df_met_agg_comb.to_excel(args.output_comb)
    print("-table saved at %s" % args.output_comb)
    print(df_met_agg_comb.to_latex(float_format="%.3f"))

# Run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draw plot showing model qualities and coefficient values.')
    parser.add_argument("--name_cohort", type=str, help="Name of the cohort.",
                        default="prism")
    parser.add_argument("--names_samples", type=str, nargs="*", help="Name of the features selection.",
                        default=[])
    parser.add_argument("--names_features", type=str, nargs="*", help="Name of the features selection.",
                        default=[])
    parser.add_argument("--names_models", type=str, nargs="*", help="Name of the model.",
                        default=[])
    parser.add_argument("--names_selections", type=str, nargs="*", help="Name of the model.",
                        default=[])
    parser.add_argument("--dir_results", type=str, help="Path to directory of results.",
                        default="../../../results/survival_analysis")
    parser.add_argument("--threshold", type=float, default=0.2,
                      help="Minimum percentage of selection of a covariate to be included inthe plot.")
    parser.add_argument("--combinations", type=str, help="Choose all or custom.", default="custom")
    parser.add_argument("--time_horizon", type=int, help="Time horizon for the c-index and brier score in days.",
                        default=180)
    parser.add_argument('--output_full', type=str, help='Path to table of scores.',
                        default="../../../results/survival_analysis/plots_prism/table_cindex_180_full.xlsx")
    parser.add_argument('--output_comb', type=str, help='Path to table of scores.',
                        default="../../../results/survival_analysis/plots_prism/table_cindex_180_comb.xlsx")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")

    main(args)
