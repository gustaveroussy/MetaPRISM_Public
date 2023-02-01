# -*- coding: utf-8 -*-
"""
@created: 28 Jun 2022
@modified: 28 Dec 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Compute numbers and percentages inside the text.
"""

import argparse
import pandas as pd
import numpy as np
import re

from pyprism.data import load_cln, load_bio, load_from_data
from pyprism.util import explode_df
from scipy.stats import mannwhitneyu, boschloo_exact
from statsmodels.sandbox.stats.multicomp import multipletests

# functions ============================================================================================================

def select_by_sample_type(df, cohort, sample_type="DNA_T\\|RNA_T", col_id="Subject_Id", cols_cln=[]):
    df_cln = load_cln(cohort)
    df_cln = df_cln.rename(columns={"Project_TCGA_More": "Tumor_Type"})

    if "DNA_P" in sample_type:
        df_cln["Sample_Type"] = df_cln["Sample_Type"].apply(lambda x: x.replace("DNA_N|DNA_T", "DNA_P"))
        mask_dnap = df_cln["Sample_Type"].apply(lambda x: "DNA_P" in x)
        cols_id_dnap = ["Sample_Id_DNA_T", "Sample_Id_DNA_N"]
        col_id_dnap = "Sample_Id_DNA_P"
        df_cln.loc[mask_dnap, col_id_dnap] = df_cln.loc[mask_dnap, cols_id_dnap].apply("_vs_".join, axis=1)

    mask = df_cln["Sample_Type"].str.contains(sample_type)
    df_cln_st = df_cln.loc[mask].copy()
    print("-INFO: selected %d/%d subjects with %s" % (sum(mask), len(mask), sample_type))

    if col_id=="Sample_Id":
        dfs_cln_sub = []
        for st in sample_type.split("\\|"):
            col_id_st = "Sample_Id_%s" % st
            df_cln_st_sub = df_cln_st.loc[~df_cln_st[col_id_st].isnull()]
            df_cln_st_sub = df_cln_st_sub.rename(columns={col_id_st: "Sample_Id"})
            dfs_cln_sub.append(df_cln_st_sub)
        df_cln_sub = pd.concat(dfs_cln_sub)
    elif col_id=="DNA_P":
        cols_sid = ["Sample_Id_DNA_T", "Sample_Id_DNA_N"]
        df_cln_sub = df_cln_st.copy()
        df_cln_sub[col_id] = df_cln_st[cols_sid].fillna("NA").apply("_vs_".join, axis=1).tolist()
    else:
        df_cln_sub = df_cln_st.copy()

    mask = df[col_id].isin(df_cln_sub[col_id].tolist())
    df_st = df.loc[mask]
    print("-INFO: selected %d/%d lines from subjects with %s" % (sum(mask), len(mask), sample_type))

    # add columns
    if len(cols_cln)>0:
        df_st = df_st.merge(df_cln_sub[[col_id]+cols_cln], how="left", on=col_id)

    return df_st


def select_by_selections(df, filepath_sam, name, col_id, col_id_sam="Sample_Id"):
    df_sam = pd.read_table(filepath_sam)
    df_sam = df_sam.loc[df_sam["Use_%s" % name]==1].copy()
    mask = df[col_id].isin(df_sam[col_id_sam])
    df = df.loc[mask].copy()
    print("-INFO: selected %d/%d lines from selected samples" % (sum(mask), len(mask)))
    return df


def numbers_abstract_intro(args):
    # percentage of patients with DNA and RNA harboring >= Tier 1 resistance
    df_alt = pd.read_table(args.CA_alt_table)
    df_cln = load_cln("prism")
    df_alt_dna_rna = select_by_sample_type(df=df_alt, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Sample_Id")
    df_cln_dna_rna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Subject_Id")

    # select samples from WES & RNAseq
    df_alt_dna_rna_sel = select_by_selections(df=df_alt_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Sample_Id", col_id_sam="Sample_Id")
    df_cln_dna_rna_sel = select_by_selections(df=df_cln_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Subject_Id", col_id_sam="Subject_Id")

    col_id = "Subject_Id"
    n_pat_tier1 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel["Res_Level_Simple"]=="Tier1"][col_id].nunique()
    n_pat_total = df_cln_dna_rna_sel[col_id].nunique()
    val = 100*n_pat_tier1/n_pat_total
    print("-VALUE: SOC resistance markers identifed in %.3g%% of patients with DNA and RNA" % val)

    n_pat_tier2 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel["Res_Level_Simple"]=="Tier2"][col_id].nunique()
    n_pat_total = df_cln_dna_rna_sel[col_id].nunique()
    val = 100*n_pat_tier2/n_pat_total
    print("-VALUE: Tier2 resistance markers identifed in %.3g%% of patients with DNA and RNA" % val)

    n_pat_tier3 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel["Res_Level_Simple"]=="Tier3"][col_id].nunique()
    n_pat_total = df_cln_dna_rna_sel[col_id].nunique()
    val = 100*n_pat_tier3/n_pat_total
    print("-VALUE: Tier3 resistance markers identifed in %.3g%% of patients with DNA and RNA" % val)

    # mask_literature = ~df_alt_dna_rna_sel["LF_Res_Drug"].isnull() & df_alt_dna_rna_sel["Res_Level_Simple"].isnull()
    # n_pat_literature = df_alt_dna_rna_sel.loc[mask_literature][col_id].nunique()
    # n_pat_total = df_cln_dna_rna_sel[col_id].nunique()
    # val = 100*n_pat_literature/n_pat_total
    # print("-VALUE: Literature resistance markers identifed in %.3g%% of patients with DNA and RNA" % val)

    # patients used for survival
    df_cln_dna_rna_sel



def numbers_part_1(args):
    df_cln = load_cln("prism")
    df_medians = pd.read_table(args.DO_violins)

    # medians survival
    print("-VALUE: medians age")
    print(df_medians[["Project_TCGA_More", "Median_Age"]].sort_values(by=["Median_Age"]))

    # tumor types and count 5 main
    col_tt = "Project_TCGA_More"
    df_cln_not_misc = df_cln.loc[df_cln[col_tt]!="MISC - Not_TCGA"]
    df_cln_not_misc_cup = df_cln.loc[~df_cln[col_tt].isin(["MISC - Not_TCGA", "Unknown_Primary"])]

    # exclude MISC - Not_TCGA from the total
    n_unique_tumor_types = df_cln_not_misc[col_tt].nunique()
    print("-VALUE: number of different tumor types %d" % n_unique_tumor_types)

    s_tt_counts = df_cln_not_misc[col_tt].value_counts()
    print("-VALUE: counts 5 main tumor types")
    print(s_tt_counts.head(5))

    # sizes of DNA and RNA cohorts

    # WES
    df_cln_dna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_T", col_id="Subject_Id")
    df_cln_not_misc_cup_dna = select_by_sample_type(df=df_cln_not_misc_cup, cohort="prism", sample_type="DNA_T",
                                                    col_id="Subject_Id")
    s_dna_tt_counts  = df_cln_not_misc_cup_dna[col_tt].value_counts()
    df_cln_dna_10 = df_cln_dna.loc[df_cln_dna[col_tt].isin(s_dna_tt_counts[s_dna_tt_counts>=10].index)]
    n_tt = df_cln_dna_10[col_tt].nunique()
    n_dna = df_cln_dna.shape[0]
    n_dna_10 = df_cln_dna_10.shape[0]
    val = 100*n_dna_10/n_dna
    print("-VALUE: DNA cohort %d types represent %.3g%% (%d/%d) of patients with DNA" % (n_tt, val, n_dna_10, n_dna))

    # RNAseq
    df_cln_rna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="RNA_T", col_id="Subject_Id")
    df_cln_not_misc_cup_rna = select_by_sample_type(df=df_cln_not_misc_cup, cohort="prism", sample_type="RNA_T",
                                                    col_id="Subject_Id")
    s_rna_tt_counts  = df_cln_not_misc_cup_rna[col_tt].value_counts()
    df_cln_rna_10 = df_cln_rna.loc[df_cln_rna[col_tt].isin(s_rna_tt_counts[s_rna_tt_counts>=10].index)]
    n_tt = df_cln_rna_10[col_tt].nunique()
    n_rna = df_cln_rna.shape[0]
    n_rna_10 = df_cln_rna_10.shape[0]
    val = 100*n_rna_10/n_rna
    print("-VALUE: RNA cohort %d types represent %.3g%% (%d/%d) of patients with RNA" % (n_tt, val, n_rna_10, n_rna))

    # WES & RNAseq
    # RNAseq
    df_cln_dna_rna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_T\\|RNA_T", col_id="Subject_Id")
    df_cln_not_misc_cup_dna_rna = select_by_sample_type(df=df_cln_not_misc_cup, cohort="prism",
                                                        sample_type="DNA_T\\|RNA_T", col_id="Subject_Id")
    s_dna_rna_tt_counts  = df_cln_not_misc_cup_dna_rna[col_tt].value_counts()
    df_cln_dna_rna_10 = df_cln_dna_rna.loc[df_cln_dna_rna[col_tt].isin(s_dna_rna_tt_counts[s_dna_rna_tt_counts>=9].index)]
    n_tt = df_cln_dna_rna_10[col_tt].nunique()
    n_dna_rna = df_cln_dna_rna.shape[0]
    n_dna_rna_10 = df_cln_dna_rna_10.shape[0]
    val = 100*n_dna_rna_10/n_dna_rna
    print("-VALUE: WES & RNA cohort %d types represent %.3g%% (%d/%d) of patients with DNA& RNA" % (n_tt, val, n_dna_rna_10, n_dna_rna))

    # metastates
    col_id = "Subject_Id"
    col_met = "Metastatic_Sites"
    df_cln_nna = df_cln.dropna(subset=[col_met])
    df_cln_met = explode_df(df_cln_nna, cols=[col_met], sep="|")

    n_pat = df_cln_nna.shape[0]
    median_met = df_cln_met.groupby([col_id]).size().median()
    n_unique_met = df_cln_met[col_met].nunique()
    print("-VALUE: median of %d metastases at %d different sites from %d patients" % (median_met, n_unique_met, n_pat))

    # biopsy sites
    df_bio = load_bio("prism")
    df_cln["Biopsy_Id"] = df_cln["Biopsy_Selected"].apply(lambda x: x.split("|")[0])
    df_cln = df_cln.merge(df_bio[["Biopsy_Id", "Biopsy_Site"]].drop_duplicates(), how="left", on="Biopsy_Id")
    s_bio_counts = df_cln["Biopsy_Site"].value_counts()
    print("-VALUE: most frequently biopsied sites")
    print(s_bio_counts.head(3)/df_cln.shape[0])

    # medians treatments
    print("-VALUE: medians treatments")
    print(df_medians[["Project_TCGA_More", "Median_Count_Drugs"]].sort_values(by=["Median_Count_Drugs"]))

    # survival fraction deceased
    df_cln_nna = df_cln.loc[~df_cln["Survival_Time"].isnull()]
    n_pat = df_cln_nna.shape[0]
    n_dec = df_cln_nna.loc[df_cln_nna["Survival_Status"]=="Deceased"].shape[0]
    print("-VALUE: %d/%d patients have died" % (n_dec, n_pat))

    # medians survival
    print("-VALUE: medians survival")
    print(df_medians[["Project_TCGA_More", "Median_Survival_Time"]].sort_values(by=["Median_Survival_Time"]))


def numbers_part_2a(args):
    df_cln = load_cln("prism")
    df_cln_dna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_P", col_id="Subject_Id")

    # mutational signatures
    df_sig = pd.read_table(args.MS_sig_table).set_index("Signature")
    df_sig = df_sig.T
    df_sig = df_sig.reset_index()
    df_sig = df_sig.rename(columns={"index": "Sample_Id_DNA_T"})
    df_cln_sig_dna = df_cln_dna.merge(df_sig, how="left", on="Sample_Id_DNA_T")

    # platinum drugs and signatures
    drugs_plat = ["CARBOPLATINE", "CISPLATINE", "OXALIPLATINE"]
    sigs_plat = ["SBS31", "SBS35"]
    col_drug = "Drugs_Before_Biopsy"

    df_cln_nna = df_cln.loc[~df_cln[col_drug].isnull()]
    df_cln_nna_plat = df_cln_nna.loc[df_cln_nna[col_drug].apply(lambda x: any(drug in x for drug in drugs_plat))]
    n_pat_nna = df_cln_nna.shape[0]
    n_pat_plat = df_cln_nna_plat.shape[0]
    print("-VALUE: %d/%d patients with drug data received platinum" % (n_pat_plat, n_pat_nna))

    df_cln_sig_dna_drug = df_cln_sig_dna.loc[~df_cln_sig_dna[col_drug].isnull()]
    df_cln_sig_dna_drug = df_cln_sig_dna_drug.loc[~df_cln_sig_dna["SBS31"].isnull()]

    for drug in drugs_plat:
        drugs_other = [x for x in drugs_plat if x!=drug]
        mask_drug = df_cln_sig_dna_drug[col_drug].apply(lambda x: drug in x)
        masks_other = [df_cln_sig_dna_drug[col_drug].apply(lambda x: drug_other in x) for drug_other in drugs_other]
        mask = mask_drug & (~masks_other[0]) & (~masks_other[1])
        print("-VALUE: number of patients that received %s only: %d" % (drug, sum(mask)))
        for sig in sigs_plat:
            mask_sig = df_cln_sig_dna_drug[sig] > 0
            prop_sig = 100*sum(mask & mask_sig) / sum(mask)
            print("-VALUE: proportion of patients with detectable sig %s: %.3g%%" % (sig, prop_sig))


    df_sum_p = pd.read_table(args.SC_summary_prism)
    df_sum_m = pd.read_table(args.SC_summary_met500)
    df_sum_t = pd.read_table(args.SC_summary_tcga)

    # add sample pair id for selection
    col_id = "Sample_Id"
    cols_sbc = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    df_sum_p[col_id] = df_sum_p[cols_sbc].fillna("NA").apply("_vs_".join, axis=1)
    df_sum_m[col_id] = df_sum_m[cols_sbc].fillna("NA").apply("_vs_".join, axis=1)
    df_sum_t[col_id] = df_sum_t[cols_sbc].fillna("NA").apply("_vs_".join, axis=1)

    # select in-design pairs
    df_sum_p = select_by_sample_type(df=df_sum_p, cohort="prism", sample_type="DNA_P", col_id="Sample_Id",
                                     cols_cln=["Tumor_Type"])
    df_sum_m = select_by_sample_type(df=df_sum_m, cohort="met500", sample_type="DNA_P", col_id="Sample_Id",
                                     cols_cln=["Tumor_Type"])
    df_sum_t = select_by_sample_type(df=df_sum_t, cohort="tcga", sample_type="DNA_P", col_id="Sample_Id",
                                     cols_cln=["Tumor_Type"])

    # select
    df_sum_p_sel = select_by_selections(df=df_sum_p, filepath_sam=args.SC_selection_prism, name="selections",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")
    df_sum_m_sel = select_by_selections(df=df_sum_m, filepath_sam=args.SC_selection_met500, name="selections",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")
    df_sum_t_sel = select_by_selections(df=df_sum_t, filepath_sam=args.SC_selection_tcga, name="selections",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")

    # compute proportion with WGD
    n_wgd_p = df_sum_p_sel.loc[df_sum_p_sel["WGD"]==1].shape[0]
    n_tot_p = df_sum_p_sel.shape[0]
    prop_p = 100*n_wgd_p/n_tot_p

    n_wgd_m = df_sum_m_sel.loc[df_sum_m_sel["WGD"]==1].shape[0]
    n_mot_m = df_sum_m_sel.shape[0]
    prop_m = 100*n_wgd_m/n_mot_m

    n_wgd_t = df_sum_t_sel.loc[df_sum_t_sel["WGD"]==1].shape[0]
    n_tot_t = df_sum_t_sel.shape[0]
    prop_t = 100*n_wgd_t/n_tot_t
    print("-VALUE: fractions of samples with WGD in PRISM: %.3g, MET500: %.3g, TCGA: %.3g" % (prop_p, prop_m, prop_t))

    # fractions of WGD by tumor type
    s_counts_wgd_p = df_sum_p_sel.groupby(["Tumor_Type"])["WGD"].mean()*100
    print("-VALUE: fractions of samples with WGD in PRISM per tumor type")
    print(s_counts_wgd_p.sort_values())

    # compare fractions using Fisher-Boschloo exact tests
    df_tests = pd.DataFrame(columns=["Cohort", "Tumor_Type", "P_Value"])

    for tt in df_sum_p_sel["Tumor_Type"].unique():
        df_sum_p_tt = df_sum_p_sel.loc[df_sum_p_sel["Tumor_Type"]==tt]
        count_wgd_tt_p = df_sum_p_tt["WGD"].sum()
        margin_wgd_tt_p = df_sum_p_tt.shape[0]
        df_sum_m_tt = df_sum_m_sel.loc[df_sum_m_sel["Tumor_Type"]==tt]
        count_wgd_tt_m = df_sum_m_tt["WGD"].sum()
        margin_wgd_tt_m = df_sum_m_tt.shape[0]

        df_sum_t_tt = df_sum_t_sel.loc[df_sum_t_sel["Tumor_Type"]==tt]
        count_wgd_tt_t = df_sum_t_tt["WGD"].sum()
        margin_wgd_tt_t = df_sum_t_tt.shape[0]

        table_p = np.array([[count_wgd_tt_p, count_wgd_tt_t],
                            [margin_wgd_tt_p-count_wgd_tt_p, margin_wgd_tt_t-count_wgd_tt_t]])
        out_test_p = boschloo_exact(table=table_p, alternative="two-sided", n=100)
        pvalue_p = out_test_p.pvalue

        table_m = np.array([[count_wgd_tt_m, count_wgd_tt_t],
                            [margin_wgd_tt_m-count_wgd_tt_m, margin_wgd_tt_t-count_wgd_tt_t]])
        out_test_m = boschloo_exact(table=table_m, alternative="two-sided", n=100)
        pvalue_m = out_test_m.pvalue

        df_tests = df_tests.append({"Cohort": "prism", "Tumor_Type": tt, "P_Value": pvalue_p}, ignore_index=True)
        df_tests = df_tests.append({"Cohort": "met500", "Tumor_Type": tt, "P_Value": pvalue_m}, ignore_index=True)

    # adjust pvalues
    df_tests["P_Value_Adj"] = multipletests(pvals=df_tests["P_Value"], alpha=0.05, method="fdr_bh")[1]

    # fractions of WGD increase PRAD
    for tt in ["PRAD", "COAD", "PAAD"]:
        df_sum_p_tt = df_sum_p_sel.loc[df_sum_p_sel["Tumor_Type"]==tt]
        prop_wgd_tt_p = 100*df_sum_p_tt["WGD"].mean()
        df_sum_m_tt = df_sum_m_sel.loc[df_sum_m_sel["Tumor_Type"]==tt]
        prop_wgd_tt_m = 100*df_sum_m_tt["WGD"].mean()
        df_sum_t_tt = df_sum_t_sel.loc[df_sum_t_sel["Tumor_Type"]==tt]
        prop_wgd_tt_t = 100*df_sum_t_tt["WGD"].mean()

        fc_wgd_tt_p = prop_wgd_tt_p/prop_wgd_tt_t
        fc_wgd_tt_m = prop_wgd_tt_m/prop_wgd_tt_t
        pvalue_tt_p = df_tests.loc[(df_tests["Tumor_Type"]==tt) & (df_tests["Cohort"]=="prism")]["P_Value_Adj"].iloc[0]
        pvalue_tt_m = df_tests.loc[(df_tests["Tumor_Type"]==tt) & (df_tests["Cohort"]=="met500")]["P_Value_Adj"].iloc[0]

        print("-VALUE: WGD fraction %s in PRISM: %.3g, MET500: %.3g" % (tt, prop_wgd_tt_p, prop_wgd_tt_m))
        print("-VALUE: fold increase WGD fraction %s in PRISM: %.2g, MET500: %.2g" % (tt, fc_wgd_tt_p, fc_wgd_tt_m))
        print("-VALUE: pvalues adj %s in PRISM: %.3g, MET500: %.3g" % (tt, pvalue_tt_p, pvalue_tt_m))

    # average genome coverage in tumors without wgd
    df_sum_p_no_wgd = df_sum_p_sel.loc[df_sum_p_sel["WGD"]==0]
    df_sum_m_no_wgd = df_sum_m_sel.loc[df_sum_m_sel["WGD"]==0]
    df_sum_t_no_wgd = df_sum_t_sel.loc[df_sum_t_sel["WGD"]==0]

    for col in ["LOSS", "GAIN"]:
        avg_value_p = 100*df_sum_p_no_wgd[col].mean()
        print("-VALUE: average genome fraction covered by %s in PRISM: %.3g%%" % (col, avg_value_p))

    # average genome coverage HLG and HD
    for col in ["LOSS:Deletion", "GAIN:HL_amplification"]:
        avg_value_p = 100*df_sum_p[col].mean()
        print("-VALUE: average genome fraction covered by %s%% in PRISM: %.3g%%" % (col, avg_value_p))

    # average number of cancer driver genes deleted/amplified
    df_cna = pd.read_table(args.SC_calls_prism)
    df_cna_ann = pd.read_table(args.SC_calls_ann_prism)
    cols_sb = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    df_cna["DNA_P"] = df_cna[cols_sb].fillna("NA").apply("_vs_".join, axis=1)
    df_cna_ann["DNA_P"] = df_cna_ann[cols_sb].fillna("NA").apply("_vs_".join, axis=1)
    df_cln = load_cln("prism")
    df_cna_sel = select_by_sample_type(df=df_cna, cohort="prism", sample_type="DNA_T", col_id="DNA_P",
                                       cols_cln=["Subject_Id", "Tumor_Type"])
    df_cna_ann_sel = select_by_sample_type(df=df_cna_ann, cohort="prism", sample_type="DNA_T", col_id="DNA_P",
                                           cols_cln=["Subject_Id", "Tumor_Type"])
    df_cln_dna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_T", col_id="Subject_Id")

    # select
    df_cna_sel = select_by_selections(df=df_cna_sel, filepath_sam=args.SC_selection_prism, name="selections",
                                      col_id="DNA_P", col_id_sam="Sample_Id")
    df_cna_ann_sel = select_by_selections(df=df_cna_ann_sel, filepath_sam=args.SC_selection_prism, name="selections",
                                          col_id="DNA_P", col_id_sam="Sample_Id")
    df_cln_dna = select_by_selections(df=df_cln_dna, filepath_sam=args.SC_selection_prism, name="selections",
                                      col_id="Subject_Id", col_id_sam="Subject_Id")


    df_cna_sel_hlg = df_cna_sel.loc[df_cna_sel["Copy_Number_More"]=="HL_gain"]
    df_cna_ann_sel_hlg = df_cna_ann_sel.loc[df_cna_ann_sel["Alteration"]=="Amplification"]
    counts_hlg = df_cna_sel_hlg.groupby("Subject_Id").size().to_frame("N_HLG").reset_index()
    counts_ann_hlg = df_cna_ann_sel_hlg.groupby("Subject_Id").size().to_frame("N_HLG_Ann").reset_index()

    df_cna_sel_hd = df_cna_sel.loc[df_cna_sel["Copy_Number_More"]=="hom_del"]
    df_cna_ann_sel_hd = df_cna_ann_sel.loc[df_cna_ann_sel["Alteration"]=="Deletion"]
    counts_hd = df_cna_sel_hd.groupby("Subject_Id").size().to_frame("N_HD").reset_index()
    counts_ann_hd = df_cna_ann_sel_hd.groupby("Subject_Id").size().to_frame("N_HD_Ann").reset_index()

    df_cln_dna = df_cln_dna.merge(counts_hlg, how="left", on="Subject_Id")
    df_cln_dna = df_cln_dna.merge(counts_ann_hlg, how="left", on="Subject_Id")
    df_cln_dna = df_cln_dna.merge(counts_hd, how="left", on="Subject_Id")
    df_cln_dna = df_cln_dna.merge(counts_ann_hd, how="left", on="Subject_Id")
    df_cln_dna["N_HLG"] = df_cln_dna["N_HLG"].fillna(0)
    df_cln_dna["N_HLG_Ann"] = df_cln_dna["N_HLG"].fillna(0)
    df_cln_dna["N_HD"] = df_cln_dna["N_HD"].fillna(0)
    df_cln_dna["N_HD_Ann"] = df_cln_dna["N_HD"].fillna(0)

    for col in ["N_HLG", "N_HD"]:
        mean_val = df_cln_dna[col].mean()
        print("-VALUE: average %s in PRISM: %.3g" % (col, mean_val))

    # check fraction of genome covered by HLG and HD
    col_id = "Subject_Id"
    cols = [col_id, "svlen"]
    length_genome = 3137144693
    df_cna_sel_hlg_len = df_cna_sel_hlg[cols].drop_duplicates().groupby(cols[:-1])[cols[-1]].sum()
    df_cna_sel_hlg_len = (df_cna_sel_hlg_len/length_genome).reset_index()
    df_cna_sel_hlg_len = df_cna_sel_hlg_len.rename(columns={"svlen": "HLG_Length"})
    df_cln_dna = df_cln_dna.merge(df_cna_sel_hlg_len, how="left", on=col_id)
    df_cln_dna["HLG_Length"] = df_cln_dna["HLG_Length"].fillna(0)

    df_cna_sel_hd_len = df_cna_sel_hd[cols].drop_duplicates().groupby(cols[:-1])[cols[-1]].sum()
    df_cna_sel_hd_len = (df_cna_sel_hd_len/length_genome).reset_index()
    df_cna_sel_hd_len = df_cna_sel_hd_len.rename(columns={"svlen": "HD_Length"})
    df_cln_dna = df_cln_dna.merge(df_cna_sel_hd_len, how="left", on=col_id)
    df_cln_dna["HD_Length"] = df_cln_dna["HD_Length"].fillna(0)

    for col in ["HLG_Length", "HD_Length"]:
        mean_val = 100*df_cln_dna[col].mean()
        print("-VALUE: average %s in PRISM: %.3g%%" % (col, mean_val))

    # get list of most frequently amplified genes among all and annotated genes
    n_pat_dna = df_cln_dna.shape[0]

    prop_hlg = 100*df_cna_sel_hlg[["Hugo_Symbol", "Chromosome"]].value_counts().to_frame("Prop")/n_pat_dna
    prop_hlg = prop_hlg.reset_index()
    prop_hlg_tt = df_cna_sel_hlg.groupby(["Hugo_Symbol", "Tumor_Type"]).size().unstack(level=-1).reset_index()
    prop_hlg = prop_hlg.merge(prop_hlg_tt, how="left", on="Hugo_Symbol")
    print(prop_hlg.head(50))
    print(prop_hlg.head(100))
    prop_hlg.loc[prop_hlg["Chromosome"]=="19"]

    prop_ann_hlg = 100*df_cna_ann_sel_hlg[["Hugo_Symbol"]].value_counts().to_frame("Prop")/n_pat_dna
    prop_ann_hlg = prop_ann_hlg.reset_index()
    prop_ann_hlg_tt = df_cna_ann_sel_hlg.groupby(["Hugo_Symbol", "Tumor_Type"]).size().unstack(level=-1).reset_index()
    prop_ann_hlg = prop_ann_hlg.merge(prop_ann_hlg_tt, how="left", on="Hugo_Symbol")
    print(prop_ann_hlg.head(10))

    # get list of most frequently deleted genes among all and annotated genes
    prop_hd = 100*df_cna_sel_hd[["Hugo_Symbol", "Chromosome"]].value_counts().to_frame("Prop")/n_pat_dna
    prop_hd = prop_hd.reset_index()
    prop_hd_tt = df_cna_sel_hd.groupby(["Hugo_Symbol", "Tumor_Type"]).size().unstack(level=-1).reset_index()
    prop_hd = prop_hd.merge(prop_hd_tt, how="left", on="Hugo_Symbol")
    print(prop_hd.head(30))

    prop_ann_hd = 100*df_cna_ann_sel_hd[["Hugo_Symbol"]].value_counts().to_frame("Prop")/n_pat_dna
    prop_ann_hd = prop_ann_hd.reset_index()
    prop_ann_hd_tt = df_cna_sel_hd.groupby(["Hugo_Symbol", "Tumor_Type"]).size().unstack(level=-1).reset_index()
    prop_ann_hd = prop_ann_hd.merge(prop_ann_hd_tt, how="left", on="Hugo_Symbol")
    print(prop_ann_hd.head(10))

    # proportion of META-PRISM WES tumors with >= 1 DNA driver event
    df_alt = pd.read_table(args.CA_alt_table)
    df_alt_dna = select_by_sample_type(df=df_alt, cohort="prism", sample_type="DNA_P", col_id="Sample_Id")
    df_alt_dna_sel = select_by_selections(df=df_alt_dna, filepath_sam=args.CA_selection, name="heatmap_dna",
                                          col_id="Sample_Id", col_id_sam="Sample_Id")
    df_cln_dna_sel = select_by_selections(df=df_cln_dna, filepath_sam=args.CA_selection, name="heatmap_dna",
                                          col_id="Subject_Id", col_id_sam="Subject_Id")

    alt_cat_mut_cna = ["Mut", "Del", "Ins", "Deletion", "Amplification"]
    df_alt_dna_sel_mut_cna = df_alt_dna_sel.loc[df_alt_dna_sel["Alteration_Category"].isin(alt_cat_mut_cna)]

    col_id = "Subject_Id"
    n_pat_alt_dna_mut_cna = df_alt_dna_sel_mut_cna[col_id].nunique()
    n_pat_cln_dna = df_cln_dna_sel[col_id].nunique()
    prop = 100*n_pat_alt_dna_mut_cna/n_pat_cln_dna
    print("-VALUE: proportion of patients in WES subcohort with >=1 DNA driver: %.3g%%" % prop)

    # frequency for top genes
    col_gen = "Hugo_Symbol"
    prop_genes = df_alt_dna_sel_mut_cna[[col_id, col_gen]].drop_duplicates()[col_gen].value_counts()/n_pat_cln_dna
    print("-VALUE: frequency driver mutation top genes")
    print(prop_genes.head(10))

    # multihit counts
    df_multi_tsg = pd.read_table(args.CA_multihit_tsg)
    df_multi_og = pd.read_table(args.CA_multihit_og)

    col_cat = "Alteration_Category"
    cat_tsg = [x for x in df_multi_tsg[col_cat] if x not in ["Mutation"]]
    print("-VALUE: fraction biallelic inactivation main TSG")
    print(df_multi_tsg.set_index(col_cat).loc[cat_tsg].sum(axis=0))
    cat_og = [x for x in df_multi_og[col_cat] if x not in ["Mutation", "Amplification"]]
    print("-VALUE: fraction multihit main OG")
    print(df_multi_og.set_index(col_cat).loc[cat_og].sum(axis=0))


def add_n_event(df_fus, df_cln, event="driver"):
    annots_known = ["ONE_PARTNER_IS_DRIVER", "Chitars_All_Human_v5.0|ONE_PARTNER_IS_DRIVER"]
    if event=="known":
        df_fus_sel = df_fus.loc[~df_fus["Annotations_Custom"].isin(annots_known)].copy()
    elif event=="driver":
        # fusions_rmv = ["TMPRSS2--ERG", "MYB--NFIB", "CLTC--VMP1", "EML4--ALK", "KIF5B--MIPOL1", "FGFR3--TACC3",
        #                "DNAJB1--PRKACA"]
        # df_fus_sel = df_fus.loc[~df_fus["Fusion_Id"].isin(fusions_rmv)].copy()
        df_fus_sel = df_fus.loc[df_fus["Annotations_Custom"].isin(annots_known)].copy()

    df_cnt_sel = df_fus_sel[["Sample_Id", "Fusion_Id"]].drop_duplicates().groupby(["Sample_Id"]).size()
    df_cnt_sel = df_cnt_sel.to_frame("N_Fus_%s" % event).reset_index()
    df_cnt_sel = df_cnt_sel.rename(columns={"Sample_Id": "Sample_Id_RNA_T"})
    df_cln_sel = df_cln.merge(df_cnt_sel, how="left", on="Sample_Id_RNA_T")
    df_cln_sel["N_Fus_%s" % event] = df_cln_sel["N_Fus_%s" % event].fillna(0).astype(int)

    return df_cln_sel


def numbers_part_2b(args):
    sheet_names = pd.ExcelFile(args.IA_heatmap).sheet_names
    dfs_ia = {sheet_name: pd.read_excel(args.IA_heatmap, sheet_name=sheet_name) for sheet_name in sheet_names}

    # percent D + F in each cohort
    cohorts = ["prism", "met500"]
    tt2tmes = {"PRAD": "D", "BLCA": "F"}
    tts = ["BLCA", "PRAD"]
    for cohort in cohorts:
        df_coh = dfs_ia["%s_percent" % cohort].set_index("Label")
        df_tcga = dfs_ia["tcga_percent"].set_index("Label")
        df_pval =  dfs_ia["%s_vs_tcga_pvals" % cohort].set_index("Label")
        for tt, tme in tt2tmes.items():
            percent_coh = df_coh.loc[tme, tt].sum()
            percent_tcga = df_tcga.loc[tme, tt].sum()
            pval = df_pval.loc[tme, tt].sum()
            print("-VALUE: fraction of TME %s in tumor type %s of %s: %.2g vs. tcga: %.2g (pval %.2g)" % \
                  (tme, tt, cohort, percent_coh, percent_tcga, pval))

    # fusions total counts
    df_cln_p = load_cln("prism").rename(columns={"Project_TCGA_More": "Tumor_Type"})
    df_cln_p_rna = select_by_sample_type(df=df_cln_p, cohort="prism", sample_type="RNA_T", col_id="Subject_Id")

    # select samples from RNAseq + fusions
    df_cln_p_fus_bur = select_by_selections(df=df_cln_p_rna, filepath_sam=args.FA_selection_prism, name="burden",
                                            col_id="Sample_Id_RNA_T", col_id_sam="Sample_Id")

    # print number of samples with fusions
    n_pat_fus_bur = df_cln_p_fus_bur.shape[0]
    n_pat_fus_can = df_cln_p_rna.loc[df_cln_p_rna["Tumor_Type"].isin(df_cln_p_fus_bur["Tumor_Type"])].shape[0]
    print("-VALUE: %d/%d RNA samples from RNAseq subcohort were analyzed for fusions" % (n_pat_fus_bur,
                                                                                         n_pat_fus_can))

    # load fusions
    df_fus_p = pd.read_table(args.FA_table_prism)
    df_fus_m = pd.read_table(args.FA_table_met500)
    df_fus_t = pd.read_table(args.FA_table_tcga)

    df_fus_p_bur = select_by_selections(df=df_fus_p, filepath_sam=args.FA_selection_prism, name="burden",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")
    df_fus_m_bur = select_by_selections(df=df_fus_m, filepath_sam=args.FA_selection_met500, name="burden",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")
    df_fus_t_bur = select_by_selections(df=df_fus_t, filepath_sam=args.FA_selection_tcga, name="burden",
                                        col_id="Sample_Id", col_id_sam="Sample_Id")

    # add number of known oncogenic fusions
    df_cln_p_fus_bur = add_n_event(df_fus=df_fus_p_bur, df_cln=df_cln_p_fus_bur, event="known")
    n_fus_known = df_cln_p_fus_bur["N_Fus_known"].sum()
    prop_fus_known  = 100*(df_cln_p_fus_bur["N_Fus_known"] > 0).mean()
    print("-VALUE: %d known oncogenic fusions in %.2g%% of PRISM RNAseq" % (n_fus_known, prop_fus_known))

    # add number of "novel" driver fusions
    df_cln_p_fus_bur = add_n_event(df_fus=df_fus_p_bur, df_cln=df_cln_p_fus_bur, event="driver")
    n_fus_driver = df_cln_p_fus_bur["N_Fus_driver"].sum()
    prop_fus_driver  = 100*(df_cln_p_fus_bur["N_Fus_driver"] > 0).mean()
    print("-VALUE: %d driver oncogenic fusions in %.2g%% of PRISM RNAseq" % (n_fus_driver, prop_fus_driver))

    # percent in top genes involved in fusions
    sheet_names = pd.ExcelFile(args.FA_heatmap_genes).sheet_names
    dfs_fa = {sheet_name: pd.read_excel(args.FA_heatmap_genes, sheet_name=sheet_name) for sheet_name in sheet_names}
    df_percent = dfs_fa["prism_count_row"].sort_values(by="Percent", ascending=False)
    df_percent["Percent"] = df_percent["Percent"]*100
    print("-VALUE: frequency fusions driver genes in PRISM RNAseq")
    print(df_percent.head(20))

    # breakpoints recurrences
    df_fus_p_bur["N_Algos_Wt_Breakpoint"] = df_fus_p_bur["Algo_Wt_Breakpoint"].apply(lambda x: len(x.split("|")))
    df_fus_p_bur["N_Algos_Wo_Breakpoint"] = df_fus_p_bur["Algo_Wo_Breakpoint"].apply(lambda x: len(x.split("|")))
    df_fus_p_bur["AR_Confidence"] = pd.Categorical(df_fus_p_bur["AR_Confidence"], categories=["medium", "high"])
    cols_order = ["N_Algos_Wt_Breakpoint", "N_Algos_Wo_Breakpoint", "AR_Confidence", "ES_Score", "PZ_Total_Reads",
                  "SF_Total_Reads"]

    genes = ["PTEN", "TP53", "RB1", "NF1"]
    dfs_fus_p_bur_gen = []
    for gene in genes:
        df_fus_p_bur_gen = df_fus_p_bur.loc[df_fus_p_bur["Fusion_Id"].apply(lambda x: gene in x)]

        # select 1 breakpoint per sample and fusion_id
        cols_keep = ["Sample_Id", "Fusion_Id", "Gene", "Breakpoint", "Direction"]
        df_fus_p_bur_gen = df_fus_p_bur_gen.sort_values(by=cols_order, ascending=False)
        df_fus_p_bur_gen_1 = df_fus_p_bur_gen.loc[df_fus_p_bur_gen["Gene_1"]==gene]
        df_fus_p_bur_gen_2 = df_fus_p_bur_gen.loc[df_fus_p_bur_gen["Gene_2"]==gene]

        df_fus_p_bur_gen_1 = df_fus_p_bur_gen_1.rename(columns={"Gene_1": "Gene", "Breakpoint_1": "Breakpoint"})
        df_fus_p_bur_gen_1["Direction"] = "1"
        df_fus_p_bur_gen_2 = df_fus_p_bur_gen_2.rename(columns={"Gene_2": "Gene", "Breakpoint_2": "Breakpoint"})
        df_fus_p_bur_gen_2["Direction"] = "2"

        df_fus_p_bur_gen_1 = df_fus_p_bur_gen_1[cols_keep].drop_duplicates(subset=cols_keep[:3], keep="first")
        df_fus_p_bur_gen_2 = df_fus_p_bur_gen_2[cols_keep].drop_duplicates(subset=cols_keep[:3], keep="first")
        df_fus_p_bur_gen = pd.concat((df_fus_p_bur_gen_1, df_fus_p_bur_gen_2), axis=0)
        df_fus_p_bur_gen["Breakpoint"] = df_fus_p_bur_gen["Breakpoint"].astype(int)
        df_fus_p_bur_gen = df_fus_p_bur_gen.sort_values(by=["Direction", "Breakpoint"])

        print("-INFO: %d fusions involving %s for %d unique samples" % (df_fus_p_bur_gen.shape[0], gene,
                                                                        df_fus_p_bur_gen["Sample_Id"].nunique()))
        dfs_fus_p_bur_gen.append(df_fus_p_bur_gen)

    df_fus_p_bur_gen = pd.concat((dfs_fus_p_bur_gen))

    # count number of recurrent breakpoint
    # a breakpoint is recurrent if it was seen in at least 2 samples, regardless of the direction
    df_fus_p_bur_gen["Breakpoint_Id"] = df_fus_p_bur_gen[["Gene", "Breakpoint"]].astype(str).apply("_".join, axis=1)
    n_breakpoints = df_fus_p_bur_gen.shape[0]
    fus_p_bur_gen_recurrent = df_fus_p_bur_gen["Breakpoint_Id"].value_counts()
    n_breakpoints_recurrent = fus_p_bur_gen_recurrent[fus_p_bur_gen_recurrent>=2].sum()
    prop_recurrent = 100*n_breakpoints_recurrent/n_breakpoints
    print("-VALUE: proportion of recurrent breakpoints for %s in PRISM: %.3g%%" % (genes, prop_recurrent))


def numbers_part_3a(args):
    # percentage of patients with DNA and RNA harboring >= Tier 1 resistance
    df_alt = pd.read_table(args.CA_alt_table)
    df_cln = load_cln("prism")
    df_alt_dna_rna = select_by_sample_type(df=df_alt, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Sample_Id")
    df_cln_dna_rna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Subject_Id")

    # select samples from WES & RNAseq
    df_alt_dna_rna_sel = select_by_selections(df=df_alt_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Sample_Id", col_id_sam="Sample_Id")
    df_cln_dna_rna_sel = select_by_selections(df=df_cln_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Subject_Id", col_id_sam="Subject_Id")

    col_id = "Subject_Id"
    for direction in ["resistance", "sensitivity"]:
        col_level = "%s_Level_Simple" % direction[:3].title()
        n_pat_tier1 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel[col_level]=="Tier1"][col_id].nunique()
        n_pat_tier2 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel[col_level]=="Tier2"][col_id].nunique()
        n_pat_tier3 = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel[col_level]=="Tier3"][col_id].nunique()
        tiers = ["Tier1", "Tier2", "Tier3"]
        n_pat_any = df_alt_dna_rna_sel.loc[df_alt_dna_rna_sel[col_level].isin(tiers)][col_id].nunique()
        n_pat_total = df_cln_dna_rna_sel[col_id].nunique()
        val_tier1 = 100*n_pat_tier1/n_pat_total
        val_any = 100*n_pat_any/n_pat_total
        print("-VALUE: SOC %s markers identifed in %.3g%% of patients with DNA and RNA" % (direction, val_tier1))
        print("-VALUE: ANY %s markers identifed in %.3g%% of patients with DNA and RNA" % (direction, val_any))

    # fraction of hypermutated tumors
    df_sum_mut = pd.read_table(args.SM_summary_mut)
    df_bed = pd.read_table(args.res_target_bed, header=None)
    df_bed["Size"] = df_bed[2] - df_bed[1] + 1
    target_size = sum(df_bed.Size)/1e6
    df_sum_mut["Burden"] = df_sum_mut["N_Mutations"]/target_size

    cols_tsb = ["Tumor_Sample_Id", "Normal_Sample_Id"]
    df_sum_mut["DNA_P"] = df_sum_mut[cols_tsb].fillna("NA").apply("_vs_".join, axis=1)
    df_sum_mut_dna_rna = select_by_sample_type(df=df_sum_mut, cohort="prism", sample_type="DNA_P\\|RNA_T",
                                               col_id="DNA_P")
    df_sum_mut_dna_rna_sel = select_by_selections(df=df_sum_mut_dna_rna, filepath_sam=args.CA_selection,
                                                  name="heatmap_all", col_id="DNA_P",
                                                  col_id_sam="Sample_Id")

    n_pat_hyper = sum(df_sum_mut_dna_rna_sel["Burden"] >= 10)
    n_pat = df_sum_mut_dna_rna_sel.shape[0]
    prop_hyper = 100*n_pat_hyper/n_pat
    print("-VALUE: fraction of hypermutation tumors %.3g%% (%d/%d) in PRISM" % (prop_hyper, n_pat_hyper, n_pat))

    # fraction of microsatellite unstable tumors
    df_msi = pd.read_table(args.MSI_prism)
    cols_tsb = ["Tumor_Sample_Id", "Normal_Sample_Id"]
    df_msi["DNA_P"] = df_msi[cols_tsb].fillna("NA").apply("_vs_".join, axis=1)
    df_msi_dna_rna = select_by_sample_type(df=df_msi, cohort="prism", sample_type="DNA_P\\|RNA_T",
                                           col_id="DNA_P")
    df_msi_dna_rna_sel = select_by_selections(df=df_msi_dna_rna, filepath_sam=args.CA_selection,
                                              name="heatmap_all", col_id="DNA_P",
                                              col_id_sam="Sample_Id")

    n_pat_msi = sum(df_msi_dna_rna_sel["Status"] == "Unstable")
    n_pat = df_msi_dna_rna_sel.shape[0]
    prop_msi = 100*n_pat_msi/n_pat
    print("-VALUE: fraction of MSI tumors %.3g%% (%d/%d) in PRISM" % (prop_msi, n_pat_msi, n_pat))

    # treatment resistances
    df_trt = pd.read_excel(args.TR_table)

    # select only patients with non-missing treatment data
    df_cln_dna_rna_sel_nna = df_cln_dna_rna_sel.loc[~df_cln_dna_rna_sel["Drugs_Before_Biopsy"].isnull()].copy()
    df_cln_dna_rna_sel_nna_exp = explode_df(df_cln_dna_rna_sel_nna, cols=["Drugs_Before_Biopsy"], sep="|")
    counts_tt = df_cln_dna_rna_sel_nna["Project_TCGA_More"].value_counts()
    df_trt = df_trt.loc[df_trt["Subject_Id"].isin(df_cln_dna_rna_sel_nna["Subject_Id"])]

    # SOC resistances
    tumor_types = df_cln_dna_rna_sel["Project_TCGA_More"].unique().tolist()
    mask_okb_tier1 = df_trt["Oncokb_Res_Level_Simple"]=="Tier1"
    mask_civ_tier1 = df_trt["Civic_Res_Level_Simple"].apply(lambda x: "Tier1" in x if type(x)==str else False)
    mask_soc = mask_okb_tier1 | mask_civ_tier1
    for tumor_type in tumor_types:
        try:
            df_dru_tt = pd.read_excel(args.TR_drugs, sheet_name=tumor_type)
            if df_dru_tt["Hide"].isnull().mean() < 1:
                df_dru_tt = df_dru_tt.loc[df_dru_tt["Hide"].str.lower()!="yes"]
            mask_class_na = df_dru_tt["Class_Plot"].isnull()
            df_dru_tt.loc[mask_class_na, "Class_Plot"] = df_dru_tt.loc[mask_class_na, "DCI"]

            mask_tt = df_trt["Tumor_Type"]==tumor_type
            mask = mask_soc & mask_tt
            if sum(mask)>0:
                df_trt_mask = df_trt.loc[mask].copy()
                trts_soc = df_trt_mask["Treatments"].unique()
                n_pat_soc = df_trt_mask["Subject_Id"].nunique()

                df_dru_tt_soc = df_dru_tt[df_dru_tt["DCI"].isin(trts_soc)].copy()
                soc_classes = df_dru_tt_soc["Class_Plot"].unique().tolist()
                for soc_class in soc_classes:
                    trts_soc_class = df_dru_tt.loc[df_dru_tt["Class_Plot"]==soc_class, "DCI"].tolist()
                    mask_cln_tt = df_cln_dna_rna_sel_nna_exp["Project_TCGA_More"]==tumor_type
                    mask_cln_trts_soc_class = df_cln_dna_rna_sel_nna_exp["Drugs_Before_Biopsy"].isin(trts_soc_class)
                    mask_cln = mask_cln_tt & mask_cln_trts_soc_class
                    n_pat_trts_soc_class = df_cln_dna_rna_sel_nna_exp.loc[mask_cln, "Subject_Id"].nunique()
                    prop_soc = 100*n_pat_soc/n_pat_trts_soc_class
                    print("-VALUE: fraction of %s tumors with SOC res %.3g%% (%d/%d) in PRISM to treatments: %s" \
                          % (tumor_type, prop_soc, n_pat_soc, n_pat_trts_soc_class, soc_class))
        except:
            pass

    # Novel resistances in LUAD EGFR INHIB 1&2 GEN
    tumor_type = "LUAD"
    for plot_class in ["EGFR INHIB. 1&2 GEN.", "EGFR INHIB. 3 GEN.", "DABRAFENIB"]:
        df_dru_tt = pd.read_excel(args.TR_drugs, sheet_name=tumor_type)
        if df_dru_tt["Hide"].isnull().mean() < 1:
            df_dru_tt = df_dru_tt.loc[df_dru_tt["Hide"].str.lower()!="yes"]
        mask_class_na = df_dru_tt["Class_Plot"].isnull()
        df_dru_tt.loc[mask_class_na, "Class_Plot"] = df_dru_tt.loc[mask_class_na, "DCI"]

        mask_dci = df_dru_tt["DCI"]==plot_class
        if sum(mask_dci) > 0:
            trts_plot_class = df_dru_tt.loc[mask_dci, "DCI"].tolist()
        else:
            trts_plot_class = df_dru_tt.loc[df_dru_tt["Class_Plot"]==plot_class, "DCI"].tolist()

        mask_cln_tt = df_cln_dna_rna_sel_nna_exp["Project_TCGA_More"]==tumor_type
        mask_cln_trts_plot_class = df_cln_dna_rna_sel_nna_exp["Drugs_Before_Biopsy"].isin(trts_plot_class)
        mask_cln = mask_cln_tt & mask_cln_trts_plot_class
        n_pat_trts_plot_class = df_cln_dna_rna_sel_nna_exp.loc[mask_cln, "Subject_Id"].nunique()

        mask_trt_trts_plot_class = df_trt["Treatments"].isin(trts_plot_class)
        mask_trt_literature = (~df_trt["LF_Gene_Alteration"].isnull()) & (df_trt["Oncokb_Res_Level_Simple"].isnull()) \
                & (df_trt["Civic_Res_Level_Simple"].isnull())
        mask_trt = mask_trt_trts_plot_class & mask_trt_literature
        n_pat_literature = df_trt.loc[mask_trt, "Subject_Id"].nunique()

        prop = 100*n_pat_literature/n_pat_trts_plot_class
        print("-VALUE: %s patients with literature-only resistance to %s: %.3g%% (%d/%d)" %
              (tumor_type, plot_class, prop, n_pat_literature, n_pat_trts_plot_class))

    # SOC fraction
    n_trt_tot  = df_trt.shape[0]
    n_trt_soc  = df_trt.loc[mask_soc].shape[0]
    prop_soc = 100*n_trt_soc/n_trt_tot
    print("-VALUE: proportion of patients-treatments with SOC resistance: %.3g%% (%d/%d)" % \
          (prop_soc, n_trt_soc, n_trt_tot))

    # Tier2 fraction
    mask_okb_tier2 = df_trt["Oncokb_Res_Level_Simple"]=="Tier2"
    mask_civ_tier2 = df_trt["Civic_Res_Level_Simple"].apply(lambda x: "Tier2" in x if type(x)==str else False)
    mask_tier2 = (mask_okb_tier2 | mask_civ_tier2) & (~mask_soc)

    n_trt_tier2  = df_trt.loc[mask_tier2].shape[0]
    prop_tier2 = 100*n_trt_tier2/n_trt_tot
    print("-VALUE: proportion of patients-treatments with Tier2 (and no better) resistance: %.3g%% (%d/%d)" % \
          (prop_tier2, n_trt_tier2, n_trt_tot))

    # Tier3 fraction
    mask_okb_tier3 = df_trt["Oncokb_Res_Level_Simple"]=="Tier3"
    mask_civ_tier3 = df_trt["Civic_Res_Level_Simple"].apply(lambda x: "Tier3" in x if type(x)==str else False)
    mask_tier3 = (mask_okb_tier3 | mask_civ_tier3) & (~mask_soc) & (~mask_tier2)

    n_trt_tier3  = df_trt.loc[mask_tier3].shape[0]
    prop_tier3 = 100*n_trt_tier3/n_trt_tot
    print("-VALUE: proportion of patients-treatments with Tier3 (and no better) resistance: %.3g%% (%d/%d)" % \
          (prop_tier3, n_trt_tier3, n_trt_tot))


    # Literature fraction
    mask_literature = (~df_trt["LF_Gene_Alteration"].isnull()) & (~mask_soc) & (~mask_tier2) & (~mask_tier3)

    n_trt_literature  = df_trt.loc[mask_literature].shape[0]
    prop_literature = 100*n_trt_literature/n_trt_tot
    print("-VALUE: proportion of patients-treatments with Literature (and no better) resistance: %.3g%% (%d/%d)" % \
          (prop_literature, n_trt_literature, n_trt_tot))

    # Add treatment category
    df_trt_cat = pd.read_excel(args.TR_drugs_cat)
    df_trt_cat = df_trt_cat.rename(columns={"DCI": "Treatments"})
    df_trt = df_trt.merge(df_trt_cat[["Treatments", "Class_Lvl_1"]].drop_duplicates(), how="left", on="Treatments")

    cat_antiangio_chemo = ["Antiangiogenic", "Chemo_Alkylating", "Chemo_Non_Alkylating"]
    prop_antiangio_chemo = df_trt["Class_Lvl_1"].isin(cat_antiangio_chemo).mean()*100
    print("-VALUE: proportion of treatments from antioangiogenic or chemo categories : %.3g%%" % prop_antiangio_chemo)


def numbers_discussion(args):
    # select PRISM WES & RNAseq cohort
    df_alt = pd.read_table(args.CA_alt_table)
    df_cln = load_cln("prism")
    df_alt_dna_rna = select_by_sample_type(df=df_alt, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Sample_Id")
    df_cln_dna_rna = select_by_sample_type(df=df_cln, cohort="prism", sample_type="DNA_P\\|RNA_T", col_id="Subject_Id")

    # select samples from WES & RNAseq
    df_alt_dna_rna_sel = select_by_selections(df=df_alt_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Sample_Id", col_id_sam="Sample_Id")
    df_cln_dna_rna_sel = select_by_selections(df=df_cln_dna_rna, filepath_sam=args.CA_selection, name="heatmap_all",
                                              col_id="Subject_Id", col_id_sam="Subject_Id")

    # percentage of alterations contributed by fusions and SCNAs
    counts_alt_cat = df_alt_dna_rna_sel["Alteration_Category"].value_counts()

    n_alt_ann = sum(counts_alt_cat)
    n_alt_scnas = counts_alt_cat["Deletion"] + counts_alt_cat["Amplification"]
    n_alt_fus = counts_alt_cat["Fusion"]
    prop_scnas = 100*n_alt_scnas/n_alt_ann
    prop_fus = 100*n_alt_fus/n_alt_ann

    print("-VALUE: proportion of annotations by SCNAs: %.3g%% (%d/%d)" % (prop_scnas, n_alt_scnas, n_alt_ann))
    print("-VALUE: proportion of annotations by fusions: %.3g%% (%d/%d)" % (prop_fus, n_alt_fus, n_alt_ann))

    # percentage of patients harboring resistances (any level) from fusions and SCNAs
    n_pat = df_cln_dna_rna_sel["Subject_Id"].nunique()
    mask_res_any = ~df_alt_dna_rna_sel["Res_Level_Simple"].isnull()
    mask_scnas = df_alt_dna_rna_sel["Alteration_Category"].isin(["Amplification", "Deletion"])
    mask_fus = df_alt_dna_rna_sel["Alteration_Category"]=="Fusion"
    n_pat_res_scnas = df_alt_dna_rna_sel.loc[mask_res_any & mask_scnas, "Subject_Id"].nunique()
    n_pat_res_fus = df_alt_dna_rna_sel.loc[mask_res_any & mask_fus, "Subject_Id"].nunique()

    prop_scnas = 100*n_pat_res_scnas/n_pat
    prop_fus = 100*n_pat_res_fus/n_pat

    print("-VALUE: proportion of patients having RES from SCNAs: %.3g%% (%d/%d)" % (prop_scnas, n_pat_res_scnas, n_pat))
    print("-VALUE: proportion of patients having RES from fusions: %.3g%% (%d/%d)" % (prop_fus, n_pat_res_fus, n_pat))


def numbers_supplementary_methods(args):
    df_bio_ind = load_bio("tcga", mode="in_design")
    df_bio_all = load_bio("tcga", mode="all")
    df_cln_ind = load_cln("tcga", mode="in_design")
    df_cln_all = load_cln("tcga", mode="all")

    mask_dna = df_cln_ind["Sample_Type"].apply(lambda x: "DNA_T" in x and "DNA_N" in x)
    df_cln_ind.loc[mask_dna, "DNA_P"] = \
            df_cln_ind.loc[mask_dna, ["Aliquot_Id_DNA_T", "Aliquot_Id_DNA_N"]].apply("_vs_".join, axis=1)

    # restrict bio tables to only selected patients
    df_bio_ind = df_bio_ind.loc[df_bio_ind["Subject_Id"].isin(df_cln_ind["Subject_Id"])]
    df_bio_all = df_bio_all.loc[df_bio_all["Subject_Id"].isin(df_cln_ind["Subject_Id"])]

    # number of patients
    print("number of patients in_design: %d" % df_cln_ind.shape[0])
    df_cln_all_not_ind = df_cln_all.loc[~df_cln_all["Subject_Id"].isin(df_cln_ind["Subject_Id"])].copy()
    df_cln_all_not_ind["Reason_Exclusion_Before"] = \
            df_cln_all_not_ind["Reason_Exclusion"].apply(lambda x: x.replace("No molecular data", ""))
    df_cln_all_not_ind["Reason_Exclusion_Before"] = \
            df_cln_all_not_ind["Reason_Exclusion_Before"].apply(lambda x: re.sub(r"\|$", "", x))
    print("reason exclusion before molecular:")
    print(df_cln_all_not_ind["Reason_Exclusion_Before"].value_counts())

    df_bio_all_dnat = df_bio_all.loc[df_bio_all["Sample_Type"]=="DNA_T"]
    mask_dnat_wxs = df_bio_all_dnat["Experimental_Strategy_Gdc"].fillna("NA").apply(lambda x: "WXS" in x) \
            | df_bio_all_dnat["Experimental_Strategy_Legacy"].fillna("NA").apply(lambda x: "WXS" in x)
    df_bio_all_dnat = df_bio_all_dnat.loc[mask_dnat_wxs]
    df_bio_all_dnan = df_bio_all.loc[df_bio_all["Sample_Type"]=="DNA_N"]
    mask_dnan_wxs = df_bio_all_dnan["Experimental_Strategy_Gdc"].fillna("NA").apply(lambda x: "WXS" in x) \
            | df_bio_all_dnan["Experimental_Strategy_Legacy"].fillna("NA").apply(lambda x: "WXS" in x)
    df_bio_all_dnan = df_bio_all_dnan.loc[mask_dnan_wxs]
    df_bio_all_dnat_m = df_bio_all_dnat[["Aliquot_Id", "Subject_Id"]].rename(columns={"Aliquot_Id": "DNA_T"})
    df_bio_all_dnan_m = df_bio_all_dnan[["Aliquot_Id", "Subject_Id"]].rename(columns={"Aliquot_Id": "DNA_N"})
    df_bio_all_dnap = df_bio_all_dnat_m.merge(df_bio_all_dnan_m, how="outer", on="Subject_Id")
    df_bio_all_dnap = df_bio_all_dnap.dropna()
    df_bio_all_dnap["DNA_P"] = df_bio_all_dnap[["DNA_T", "DNA_N"]].apply("_vs_".join, axis=1)
    df_bio_all_dnap = df_bio_all_dnap.merge(df_cln_all[["Subject_Id", "Gender"]], how="left", on="Subject_Id")

    # somatic mutations
    # get number of pairs not analyzed/excluded by MC3
    df_maf_sam = load_from_data("tcga/wes/somatic_maf/sample_list.tsv")
    df_maf_sam = df_maf_sam.rename(columns={"Tumor_Sample_Id": "DNA_T", "Normal_Sample_Id": "DNA_N"})
    cols_sam = ["DNA_T", "DNA_N", "Comment", "QC"]
    df_bio_all_dnap = df_bio_all_dnap.merge(df_maf_sam[cols_sam], how="left", on=["DNA_T", "DNA_N"])
    df_bio_all_dnap["QC"] = df_bio_all_dnap["QC"].fillna("FAIL")
    mask_fail = df_bio_all_dnap["QC"]=="FAIL"
    df_bio_all_dnap.loc[mask_fail, "Comment"] = df_bio_all_dnap.loc[mask_fail, "Comment"].fillna("MC3 not analyzed")

    mask_1 = ~df_bio_all_dnap["Comment"].apply(lambda x: "MC3" in x if type(x)==str else False)
    mask_2 = ~df_bio_all_dnap["Comment"].apply(lambda x: "GDC Contamination" in x if type(x)==str else False)
    mask_3 = ~df_bio_all_dnap["DNA_P"].isin(df_cln_ind["DNA_P"])

    print("number of pairs: %d" % df_bio_all_dnap.shape[0])
    print("mask_1: %d"  % sum(mask_1))
    print("not mask_1: %d"  % sum(~mask_1))
    print("mask_1 mask_2: %d"  % sum(mask_1 & mask_2))
    print("mask_1 not mask_2: %d"  % sum(mask_1 & ~mask_2))
    print("mask_1 mask_2 mask 3: %d"  % sum(mask_1 & mask_2 & mask_3))
    print("mask_1 mask_2 not mask_3: %d"  % sum(mask_1 & mask_2 & ~mask_3))

    # somatic CNAs
    del df_bio_all_dnap["QC"]
    del df_bio_all_dnap["Comment"]

    df_bam_avail = pd.read_csv("../../../data/tcga/clinical/raw_bio_files/tcga_wxs_bams_isb_cgc.csv")
    df_bam_avail = df_bam_avail.rename(columns={"aliquot_barcode": "Aliquot_Id"})
    df_bam_avail["BAM_Avail"] = "Yes"
    df_bam_avail["DNA_T"] = df_bam_avail["Aliquot_Id"]
    df_bam_avail["BAM_DNA_T"] = "Yes"
    df_bam_avail["DNA_N"] = df_bam_avail["Aliquot_Id"]
    df_bam_avail["BAM_DNA_N"] = "Yes"
    df_bio_all_dnap = df_bio_all_dnap.merge(df_bam_avail[["DNA_T", "BAM_DNA_T"]], on="DNA_T", how="left")
    df_bio_all_dnap = df_bio_all_dnap.merge(df_bam_avail[["DNA_N", "BAM_DNA_N"]], on="DNA_N", how="left")

    df_tnp_facets = pd.read_table("../../../data/tcga/wes/config_facets/tumor_normal_pairs.12119.tsv")
    df_tnp_facets["Run"] = "SUCCESS"
    df_bio_all_dnap = df_bio_all_dnap.merge(df_tnp_facets[["DNA_P", "Run"]].drop_duplicates(), how="left")

    df_cna_sam = load_from_data("tcga/wes/somatic_cna/sample_list.tsv")
    df_cna_sam["DNA_P"] = df_cna_sam[["Tumor_Sample_Id", "Normal_Sample_Id"]].apply("_vs_".join, axis=1)
    df_bio_all_dnap = df_bio_all_dnap.merge(df_cna_sam[["DNA_P", "Comment", "QC"]], how="left", on="DNA_P")

    mask_1 = df_bio_all_dnap[["BAM_DNA_T", "BAM_DNA_N", "Gender"]].isnull().sum(axis=1)==0
    mask_2 = df_bio_all_dnap["Run"]=="SUCCESS"
    mask_3 = ~df_bio_all_dnap["Comment"].apply(lambda x: "GDC Contamination" in x if type(x)==str else False)
    mask_4 = ~df_bio_all_dnap["Comment"].apply(lambda x: "MP Failed QC" in x if type(x)==str else False)
    mask_5 = df_bio_all_dnap["DNA_P"].isin(df_cln_ind["DNA_P"])
    print("shape: %d"  % df_bio_all_dnap.shape[0])
    print("mask_1: %d"  % sum(mask_1))
    print("not mask_1: %d"  % sum(~mask_1))
    print("mask_1 mask_2: %d"  % sum(mask_1 & mask_2))
    print("mask_1 not mask_2: %d"  % sum(mask_1 & ~mask_2))
    print("mask_1 mask_2 mask 3: %d"  % sum(mask_1 & mask_2 & mask_3))
    print("mask_1 mask_2 not mask_3: %d"  % sum(mask_1 & mask_2 & ~mask_3))
    print("mask_1 mask_2 mask 3 mask 4: %d"  % sum(mask_1 & mask_2 & mask_3 & mask_4))
    print("mask_1 mask_2 mask_3 not mask 4: %d"  % sum(mask_1 & mask_2 & mask_3 & ~mask_4))
    print("mask_1 mask_2 mask 3 mask 4 mask 5: %d"  % sum(mask_1 & mask_2 & mask_3 & mask_4 & mask_5))
    print("mask_1 mask_2 mask_3 mask 4 not mask 5: %d"  % sum(mask_1 & mask_2 & mask_3 & mask_4 & ~mask_5))

    # tumor msi
    del df_bio_all_dnap["QC"]
    del df_bio_all_dnap["Comment"]

    df_msi_sam = load_from_data("tcga/wes/somatic_msi/sample_list.tsv")
    df_msi_sam["DNA_P"] = df_msi_sam[["Tumor_Sample_Id", "Normal_Sample_Id"]].apply("_vs_".join, axis=1)
    df_msi_sam["Run_msi"] = "Yes"
    df_bio_all_dnap = df_bio_all_dnap.merge(df_msi_sam[["DNA_P", "Run_msi", "Comment", "QC"]], how="left", on="DNA_P")

    mask_1 = df_bio_all_dnap["Run_msi"]=="Yes"
    mask_2 = ~df_bio_all_dnap["Comment"].apply(lambda x: "GDC Contamination" in x if type(x)==str else False)
    mask_3 = df_bio_all_dnap["DNA_P"].isin(df_cln_ind["DNA_P"].dropna())
    print("shape: %d"  % df_bio_all_dnap.shape[0])
    print("mask_1: %d"  % sum(mask_1))
    print("not mask_1: %d"  % sum(~mask_1))
    print("mask_1 mask_2: %d"  % sum(mask_1 & mask_2))
    print("mask_1 not mask_2: %d"  % sum(mask_1 & ~mask_2))
    print("mask_1 mask_2 mask 3: %d"  % sum(mask_1 & mask_2 & mask_3))
    print("mask_1 mask_2 not mask_3: %d"  % sum(mask_1 & mask_2 & ~mask_3))

    # tumor gene expression
    df_bio_all_rnat = df_bio_all.loc[df_bio_all["Sample_Type"]=="RNA_T"]
    mask_1 = df_bio_all_rnat["RNA_kallisto-tximport"]==1
    mask_2 = df_bio_all_rnat["RNA_kallisto-tximport_qc_status"]=="PASS"
    mask_3 = df_bio_all_rnat["Aliquot_Id"].isin(df_cln_ind["Aliquot_Id_RNA_T"].dropna())
    print("shape: %d"  % df_bio_all_rnat.shape[0])
    print("mask_1: %d"  % sum(mask_1))
    print("not mask_1: %d"  % sum(~mask_1))
    print("mask_1 mask_2: %d"  % sum(mask_1 & mask_2))
    print("mask_1 not mask_2: %d"  % sum(mask_1 & ~mask_2))
    print("mask_1 mask_2 mask 3: %d"  % sum(mask_1 & mask_2 & mask_3))
    print("mask_1 mask_2 not mask_3: %d"  % sum(mask_1 & mask_2 & ~mask_3))

    # germline mutations
    df_bio_all_dnan = df_bio_all.loc[df_bio_all["Sample_Type"]=="DNA_N"]
    mask_gdc_wxs = df_bio_all_dnan["Experimental_Strategy_Gdc"].apply(lambda x: "WXS" in x if type(x)==str else False)
    mask_leg_wxs = df_bio_all_dnan["Experimental_Strategy_Legacy"].apply(lambda x: "WXS" in x if type(x)==str else False)
    df_bio_all_dnan = df_bio_all_dnan.loc[mask_gdc_wxs | mask_leg_wxs]
    mask_1 = df_bio_all_dnan["WES_germline_maf"]==1
    mask_2 = df_bio_all_dnan["WES_germline_maf_qc_comment"].isnull()
    mask_3 = df_bio_all_dnan["Aliquot_Id"].isin(df_cln_ind["Aliquot_Id_DNA_N"].dropna())
    print("shape: %d"  % df_bio_all_dnan.shape[0])
    print("mask_1: %d"  % sum(mask_1))
    print("not mask_1: %d"  % sum(~mask_1))
    print("mask_1 mask_2: %d"  % sum(mask_1 & mask_2))
    print("mask_1 not mask_2: %d"  % sum(mask_1 & ~mask_2))
    print("mask_1 mask_2 mask 3: %d"  % sum(mask_1 & mask_2 & mask_3))
    print("mask_1 mask_2 not mask_3: %d"  % sum(mask_1 & mask_2 & ~mask_3))


def main(args):
    print("="*40)
    print("ABSTRACT/INTRO")
    numbers_abstract_intro(args)
    print("="*40)

    print("="*40)
    print("PART 1")
    numbers_part_1(args)
    print("="*40)

    print("="*40)
    print("PART 2A")
    numbers_part_2a(args)
    print("="*40)

    print("="*40)
    print("PART 2B")
    numbers_part_2b(args)
    print("="*40)

    print("="*40)
    print("PART 3A")
    numbers_part_3a(args)
    print("="*40)



# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute numbers and percentages in the text.')
    parser.add_argument('--CA_alt_table', type=str, help='Table of aggregated alterations from CA.',
                        default="../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv")
                # default="../../../results/treatment_resistances/annotations/aggregated_alterations_prism_all_annot.tsv")
    parser.add_argument('--CA_selection', type=str, help='Table of sample selections from CA.',
                        default="../../../results/combined_alterations/selection/selection_samples_prism.tsv")
    parser.add_argument('--CA_multihit_og', type=str, help='Table of oncogenes multihit counts from CA.',
                        default="../../../results/combined_alterations/other_plots/table_multihits_oncogenes.tsv")
    parser.add_argument('--CA_multihit_tsg', type=str, help='Table of tumorsuppressors multihit counts from CA.',
                        default="../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors.tsv")
    parser.add_argument('--DO_violins', type=str, help='Table of violins medians from DO.',
                        default="../../../results/data_overview/violins/violins_medians.tsv")
    parser.add_argument('--IA_subtypes', type=str, help='Table of TME subtypes from IA.',
                        default="../../../results/immuno_analysis/prism/tables/mfp_subtypes_predicted_LogisticRegression.tsv")
    parser.add_argument('--IA_heatmap', type=str, help='Table underlying the heatmap of TME.',
                        default="../../../results/immuno_analysis/heatmap/tables.xlsx")
    parser.add_argument('--MS_sig_table', type=str, help='Table of mutational signatures from MS.',
                        default="../../../results/mutational_signatures/projection_known_signatures/MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_prism.tsv")
    parser.add_argument('--SC_summary_prism', type=str, help='Table of SCNA summary stats for prism from SC.',
                        default="../../../data/prism/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz")
    parser.add_argument('--SC_summary_met500', type=str, help='Table of SCNA summary stats for met500 from SC.',
                        default="../../../data/met500/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz")
    parser.add_argument('--SC_summary_tcga', type=str, help='Table of SCNA summary stats for met500 from SC.',
                        default="../../../data/tcga/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz")
    parser.add_argument('--SC_calls_prism', type=str, help='Table of SCNA calls for prism from SC.',
                        default="../../../data/prism/wes/somatic_cna/somatic_calls.tsv.gz")
    parser.add_argument('--SC_calls_ann_prism', type=str, help='Table of SCNA calls for prism from SC.',
                        default="../../../data/prism/wes/somatic_cna/somatic_calls_union_ann.tsv.gz")
    parser.add_argument('--SC_selection_prism', type=str, help='Table of samples selection for prism from SC.',
                        default="../../../results/somatic_cnas/selection/selection_samples_prism.tsv")
    parser.add_argument('--SC_selection_met500', type=str, help='Table of samples selection for met500 from SC.',
                        default="../../../results/somatic_cnas/selection/selection_samples_met500.tsv")
    parser.add_argument('--SC_selection_tcga', type=str, help='Table of samples selection for tcga from SC.',
                        default="../../../results/somatic_cnas/selection/selection_samples_tcga.tsv")
    parser.add_argument('--FA_table_prism', type=str, help='Table of fusions for prism from FA.',
                        default="../../../data/prism/rna/fusions/prism_annotated_filtered.tsv.gz")
    parser.add_argument('--FA_table_met500', type=str, help='Table of fusions for met500 from FA.',
                        default="../../../data/met500/rna/fusions/met500_annotated_filtered.tsv.gz")
    parser.add_argument('--FA_table_tcga', type=str, help='Table of fusions for tcga from FA.',
                        default="../../../data/tcga/rna/fusions/tcga_annotated_filtered.tsv.gz")
    parser.add_argument('--FA_selection_prism', type=str, help='Table of samples selection for prism from FA.',
                        default="../../../results/fusions_analysis/selection/selection_samples_prism.tsv")
    parser.add_argument('--FA_selection_met500', type=str, help='Table of samples selection for met500 from FA.',
                        default="../../../results/fusions_analysis/selection/selection_samples_met500.tsv")
    parser.add_argument('--FA_selection_tcga', type=str, help='Table of samples selection for tcga from FA.',
                        default="../../../results/fusions_analysis/selection/selection_samples_tcga.tsv")
    parser.add_argument('--FA_heatmap_genes', type=str, help='Table of samples selection for tcga from FA.',
                        default="../../../results/fusions_analysis/heatmap/tables_genes.xlsx")
    parser.add_argument('--SM_summary_mut', type=str, help='Table of mutations summary for PRISM from SM.',
                        default="../../../data/prism/wes/summary/somatic_maf.tsv")
    parser.add_argument('--MSI_prism', type=str, help='Table of MSI for PRISM.',
                        default="../../../data/prism/wes/somatic_msi/somatic_msi.tsv")
    parser.add_argument('--res_target_bed', type=str, help='BED of targets used in PRISM.',
                        default="../../../data/resources/target_files/all_targets_intersect_padded_10n.bed")
    parser.add_argument('--TR_table', type=str, help='Table of treatment resistances for PRISM from TR.',
                        default="../../../results/treatment_resistances/annotations/treatment_resistances_prism.xlsx")
    parser.add_argument('--TR_drugs', type=str, help='Table of treatment resistances classes for PRISM from TR.',
                        default="../../../data/resources/drug_tables/Table_Drugs_Groups_Oncoplot_29072022.xlsx")
    parser.add_argument('--TR_drugs_cat', type=str, help='Table of drug categories.',
                        default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
