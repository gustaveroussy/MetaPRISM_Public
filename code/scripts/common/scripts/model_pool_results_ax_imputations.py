# -*- coding: utf-8 -*-
"""
@created: Dec 04 2021
@modified: May 06 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Pool coefficient/standard error and quality metrics estimates across imputed tables using Rubin's formula.
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy.stats import mode


def load_met_cov_across_imputations(run_results, n_imputations):
    dfs_met = []
    dfs_cov = []

    for m in range(1,n_imputations+1):
        filename_met = "mets.imputed_%d.tsv.gz" % m
        filename_cov = "covs.imputed_%d.tsv.gz" % m

        filepath_met = os.path.join(run_results, filename_met)
        filepath_cov = os.path.join(run_results, filename_cov)

        dfs_met.append(pd.read_table(filepath_met))
        dfs_cov.append(pd.read_table(filepath_cov))

    return pd.concat(dfs_met), pd.concat(dfs_cov)


def load_met_cov_complete_cases(run_results):
    filename_met = "mets.complete_cases.tsv.gz"
    filename_cov = "covs.complete_cases.tsv.gz"

    filepath_met = os.path.join(run_results, filename_met)
    filepath_cov = os.path.join(run_results, filename_cov)

    return pd.read_table(filepath_met), pd.read_table(filepath_cov)


def nanmode(x):
    x_nna = x.dropna()
    if len(x_nna)==0:
        return np.nan
    else:
        return mode(x.dropna()).mode[0]


def poolstr(x):
    x_l = x.dropna().tolist()
    x_uni_all = sorted(list(set(x_l)))
    if len(x_uni_all)==0:
        return np.nan
    else:
        x_cnt_all = [x_l.count(x_uni) for x_uni in x_uni_all]
        return ";".join(["%s (x%d)" % (x_uni, x_cnt) for x_uni, x_cnt in zip(x_uni_all, x_cnt_all)])


def pool_point_estimate(df, cols_gby):
    cols_avoid = ["Table"]
    cols_str = ["Warning", "Error", "Warning_Selection"]
    cols_str = [x for x in cols_str if x in df]
    cols_est = [x for x in df if x not in cols_str + cols_gby + cols_avoid]
    dtypes = df[cols_est].dtypes
    mask_num = dtypes.apply(lambda x: np.issubdtype(x,np.number))
    cols_est_num = dtypes[mask_num].index.tolist()
    cols_est_obj = dtypes[~mask_num].index.tolist()

    dt_agg_num = {col_est: np.nanmean for col_est in cols_est_num}
    dt_agg_obj = {col_est: nanmode for col_est in cols_est_obj}
    dt_agg_str = {col_str: poolstr for col_str in cols_str}
    dt_agg = {**dt_agg_num, **dt_agg_obj, **dt_agg_str}

    return df.groupby(cols_gby).agg(dt_agg).reset_index(drop=False)


def pool_variance_estimate(df, n_imputations, col_point="Coefficient", col_std="Standard_Error",
                           cols_gby=["Covariate", "Run", "Selection", "Model", "Features"]):
    col_var = "%s_Squared" % col_std
    df[col_var] = df[col_std].apply(np.square)

    col_var_w = "%s_Within" % col_var
    col_var_b = "%s_Between" % col_var

    # within-imputation variance
    cols_est_1 = [col_var]
    dt_agg_1 = {col_est: np.nanmean for col_est in cols_est_1}
    df_wimp_var = df.groupby(cols_gby).agg(dt_agg_1).rename(columns={col_var: col_var_w})

    # between-imputation variance
    cols_est_2 = [col_point]
    dt_agg_2 = {col_est: lambda x: np.nanvar(x, ddof=1) for col_est in cols_est_2}
    df_bimp_var = df.groupby(cols_gby).agg(dt_agg_2).rename(columns={col_point: col_var_b})

    # total variance estimate
    df_pool = pd.concat((df_bimp_var, df_wimp_var), axis=1)
    df_pool[col_var] = df_pool[col_var_w] + (1+1/n_imputations)*df_pool[col_var_b]
    df_pool[col_std] = df_pool[col_var].apply(np.sqrt)
    df_pool["%s_Relative_Increase" % col_var] = (1+1/n_imputations)*df_pool[col_var_b]/df_pool[col_var_w]
    del df_pool[col_var_b]
    del df_pool[col_var_w]

    return df_pool.reset_index()


def merge_on_common_cols(df_a, df_b, **kwargs):
    cols_c = list(set(df_a.columns).intersection(set(df_b.columns)))
    return df_a.merge(df_b, on=cols_c, **kwargs)


def main(args):
    imputation_performed = any([x.startswith("mets.imputed") for x in os.listdir(args.run_results)])

    if imputation_performed:
        df_met, df_cov = load_met_cov_across_imputations(run_results=args.run_results,
                                                         n_imputations=args.n_imputations)
    else:
        df_met, df_cov = load_met_cov_complete_cases(run_results=args.run_results)

    # split cross-validation or bootstrap/jackknife estimates from full dataset estimates
    # split bootstrap/jackknife estimates from full dataset estimates
    df_cov_rept = df_cov.loc[df_cov["Run"]!="Full"].copy()
    df_cov_full = df_cov.loc[df_cov["Run"]=="Full"].copy()
    del df_cov

    # harmonize NAs across metrics
    cols_met = ["F_comp", "C_simple", "C_ipcw", "Brier",
                "accuracy", "recall", "precision", "f1_weighted", "f1_macro", "roc_auc"]
    cols_met = [x for x in cols_met if x in df_met]

    mask_nas = df_met[cols_met].isnull().sum(axis=1) > 0
    df_met.loc[mask_nas, cols_met] = np.nan

    # pool metric estimates across imputations by mean
    cols_gby_met = ["Time", "Split", "Run", "Grid_Run", "Selection", "Model", "Features"]
    cols_gby_met = [x for x in cols_gby_met if x in df_met]

    df_met[cols_gby_met] = df_met[cols_gby_met].fillna("N/A")
    df_met_pool = pool_point_estimate(df_met, cols_gby=cols_gby_met)

    # pool point estimates across imputations using Rubin's formula (mean) 
    cols_gby_cov = ["Covariate", "Run", "Grid_Run", "Selection", "Model", "Features"]
    cols_gby_cov = [x for x in cols_gby_cov if x in df_cov_rept]

    df_cov_pool_rept = pool_point_estimate(df=df_cov_rept, cols_gby=cols_gby_cov)
    df_cov_pool_full = pool_point_estimate(df=df_cov_full, cols_gby=cols_gby_cov)

    if "Standard_Error" in df_cov_rept and imputation_performed:
        df_cov_pool_rept_std = pool_variance_estimate(df_cov_rept, n_imputations=args.n_imputations,
                                                      col_point="Coefficient", col_std="Standard_Error",
                                                      cols_gby=cols_gby_cov)

        df_cov_pool_full_std = pool_variance_estimate(df_cov_full, n_imputations=args.n_imputations,
                                                      col_point="Coefficient", col_std="Standard_Error",
                                                      cols_gby=cols_gby_cov)

        del df_cov_pool_rept["Standard_Error"]
        del df_cov_pool_full["Standard_Error"]

        df_cov_pool_rept = merge_on_common_cols(df_cov_pool_rept, df_cov_pool_rept_std)
        df_cov_pool_full = merge_on_common_cols(df_cov_pool_full, df_cov_pool_full_std)

    # put back reptstrap estimates with full dataset estimates
    df_cov_pool = pd.concat((df_cov_pool_rept, df_cov_pool_full), axis=0)

    # save
    df_met_pool.replace("N/A", np.nan).to_csv(args.output_met, index=False, sep="\t")
    print("-file saved at %s" % args.output_met)
    df_cov_pool.replace("N/A", np.nan).to_csv(args.output_cov, index=False, sep="\t")
    print("-file saved at %s" % args.output_cov)


# Run ==================================================================================================================

if __name__ == "__main__":
    results_folder = "../../../results/survival_analysis"
    default_folder = "%s/models/BRCA_dna_and_rna/cln_biol_2/lasso_coxph_standard" % results_folder

    parser = argparse.ArgumentParser(description='Pool model quality metrics and fitted coefficients across imputations.')
    parser.add_argument("--run_results", type=str, help="Path to directory of run results.", default=default_folder)
    parser.add_argument("--n_imputations", type=int, default=2, help="Number of imputed tables.")
    parser.add_argument('--output_met', type=str, help='Path to output table of pooled model quality metrics.',
                        default="%s/mets.pooled_none_coxph_standard_ax_imp.tsv.gz" % default_folder)
    parser.add_argument('--output_cov', type=str, help='Path to output table of pooled coefficient estimates.',
                        default="%s/covs.pooled_none_coxph_standard_ax_imp.tsv.gz" % default_folder)
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
