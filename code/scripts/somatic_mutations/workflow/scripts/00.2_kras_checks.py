# -*- coding: utf-8 -*-
"""
@created: Mar 15 2022
@modified: Jul 06 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Analyse mutations data for PAAD from different sources and detail discrepant mutation calls for KRAS.
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys

from pyprism.data import load_cln, load_bio, load_wes_mut

# functions ============================================================================================================

def convert_to_str(x):
    try:
        if x==int(x):
            y="%d" % int(x)
        else:
            raise
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = np.nan
        except:
            y = x
    return y

def get_sam_awg(sam_awg):
    filepath_sam_awg = sam_awg
    df_sam_awg = pd.read_excel(filepath_sam_awg, sheet_name="FreezeSamples", skiprows=1)
    df_sam_awg = df_sam_awg.rename(columns={"Tumor Sample ID": "Biopsy_Vial_Id"})

    return df_sam_awg


def get_sam_bio():
    df_bio = load_bio(study="tcga", mode="all")
    df_cln = load_cln(study="tcga", mode="all")
    df_bio = df_bio.merge(df_cln[["Subject_Id", "Project_TCGA_More"]], how="left", on="Subject_Id")

    mask_paad = df_bio["Project_TCGA_More"]=="PAAD"
    mask_dnat = df_bio["Sample_Type"]=="DNA_T"
    mask_dnan = df_bio["Sample_Type"]=="DNA_N"
    mask_raws = df_bio["Data_Category_Legacy"]=="Raw sequencing data"
    mask_mc3 = df_bio["WES_somatic_maf"]==1

    cols = ["Aliquot_Id", "Sample_Id", "Biopsy_Id", "Biopsy_Vial_Id", "Analyte_Type", "Analyte_Type_Code",
            "Slide_Percent_Tumor_Cells"]
    df_sam_bio = df_bio.loc[mask_paad & mask_dnat & mask_raws & mask_mc3][cols].drop_duplicates()

    return df_sam_bio


def get_sam_alt(sam_alt):
    filepath_sam_alt = sam_alt
    df_alt = pd.read_table(filepath_sam_alt)
    df_cln = load_cln(study="tcga")
    df_bio = load_bio("tcga")
    df_bio = df_bio.merge(df_cln[["Subject_Id", "Project_TCGA_More"]], how="left", on="Subject_Id")
    df_alt = df_alt.merge(df_bio[["Sample_Id","Project_TCGA_More", "Sample_Type"]].drop_duplicates(), how='left',
                          on="Sample_Id")
    df_alt = df_alt.loc[(df_alt["Use_heatmap"]==1)&(df_alt["Sample_Type"]=="DNA_T")]

    mask_paad = df_alt["Project_TCGA_More"]=="PAAD"
    df_sam_alt = df_alt.loc[mask_paad][["Sample_Id"]].drop_duplicates()
    df_sam_alt["Biopsy_Vial_Id"] = df_sam_alt["Sample_Id"].str[:16]

    return df_sam_alt


def get_sam_mc3_pub(mc3_pub):
    filepath_sam_mc3 = mc3_pub
    df_sam_mc3 = pd.read_table(filepath_sam_mc3, usecols=["Tumor_Sample_Barcode"])
    df_sam_mc3["Biopsy_Vial_Id"] = df_sam_mc3["Tumor_Sample_Barcode"].str[:16]

    return df_sam_mc3


def get_sam_mc3_con(mc3_con):
    filepath_sam_mc3 = mc3_con
    df_sam_mc3 = pd.read_table(filepath_sam_mc3, usecols=["Tumor_Sample_Barcode"])
    df_sam_mc3["Biopsy_Vial_Id"] = df_sam_mc3["Tumor_Sample_Barcode"].str[:16]

    return df_sam_mc3


def outer_join_sam_tables(df_sam_awg, df_sam_bio, df_sam_alt, df_sam_mc3_pub, df_sam_mc3_con):
    # prepare results table
    cols_sam_awg = [col_id, "Initial Slide Tumor Cellularity", "Pathologist Reviewed Tumor Cellularity",
                    "ABSOLUTE Purity", "Ploidy", "KRAS Mutated (1 or 0)"]
    df_res = df_sam_awg[cols_sam_awg]
    df_res["In_TCGA_AWG_2017"] = "Yes"

    # add sam bio
    cols_sam_bio = [col_id, "In_TCGA_Bio", "Slide_Percent_Tumor_Cells"]
    df_sam_bio["In_TCGA_Bio"] = "Yes"
    df_res = df_res.merge(df_sam_bio[cols_sam_bio].drop_duplicates(), how="outer", on=col_id)

    # add sam alt
    cols_sam_alt = [col_id, "In_Htmp_Alt"]
    df_sam_alt["In_Htmp_Alt"] = "Yes"
    df_res = df_res.merge(df_sam_alt[cols_sam_alt].drop_duplicates(), how="outer", on=col_id)

    # add sam mc3 pub
    cols_sam_mc3_pub = [col_id, "In_MC3_Public"]
    df_sam_mc3_pub["In_MC3_Public"] = "Yes"
    df_res = df_res.merge(df_sam_mc3_pub[cols_sam_mc3_pub].drop_duplicates(), how="outer", on=col_id)

    # add sam mc3 con
    cols_sam_mc3_con = [col_id, "In_MC3_Controlled"]
    df_sam_mc3_con["In_MC3_Controlled"] = "Yes"
    df_res = df_res.merge(df_sam_mc3_con[cols_sam_mc3_con].drop_duplicates(), how="outer", on=col_id)

    return df_res


def get_mut_bio():
    df_sam_bio = get_sam_bio()
    df_mut_fil = load_wes_mut(study="tcga", mode="somatic_filters")

    mask_bio = df_mut_fil["Tumor_Sample"].isin(df_sam_bio["Aliquot_Id"])
    df_mut_fil = df_mut_fil.loc[mask_bio]
    df_mut_fil.loc[df_mut_fil["TYPE"]=="INDEL", "POS"] = df_mut_fil.loc[df_mut_fil["TYPE"]=="INDEL", "POS"] + 1
    df_mut_fil["Biopsy_Vial_Id"] = df_mut_fil["Tumor_Sample"].str[:16]

    df_mut_pas = load_wes_mut(study="tcga", mode="somatic_maf")
    mask_bio = df_mut_pas["Tumor_Sample_Barcode"].isin(df_sam_bio["Aliquot_Id"])
    df_mut_pas = df_mut_pas.loc[mask_bio]
    df_mut_pas["Biopsy_Vial_Id"] = df_mut_pas["Tumor_Sample_Barcode"].str[:16]

    return df_mut_fil, df_mut_pas


def get_mut_mc3_con(mc3_con):
    filepath_mut_mc3 = mc3_con
    df_mut_mc3 = pd.read_table(filepath_mut_mc3)
    df_mut_mc3["Biopsy_Vial_Id"] = df_mut_mc3["Tumor_Sample_Barcode"].str[:16]

    return df_mut_mc3


def add_gene_info_mut(df_res, df_mut, col_id, gene, suffix, cols_show):
    df_mut_sub = df_mut.loc[df_mut["Hugo_Symbol"]==gene].copy()

    if all([col in cols_show for col in ["t_alt_count", "t_depth", "t_vaf"]]):
        df_mut_sub["t_vaf"] = df_mut_sub["t_alt_count"]/df_mut_sub["t_depth"]

    for col_show in cols_show:
        df_mut_sub[col_show] = df_mut_sub[col_show].apply(convert_to_str).astype(str)
    df_mut_info = df_mut_sub.groupby(col_id).agg({col_show: " | ".join for col_show in cols_show}).reset_index()

    if "HGVSp_Short" in cols_show:
        df_mut_info["HGVSp_Short"] = df_mut_info["HGVSp_Short"].apply(lambda x: x.replace("%3D","=") if type(x)==str else x)

    old2new = {x: "%s_%s" % (x, suffix) for x in cols_show}
    df_mut_info = df_mut_info.rename(columns=old2new)
    df_res = df_res.merge(df_mut_info, how="outer", on=col_id)

    return df_res


def main(args):
    col_id = "Biopsy_Vial_Id"

    # supplementary PAAD paper TCGA AWG 2017
    df_sam_awg = get_sam_awg(args.sam_awg)

    # biospecimen all for TCGA
    df_sam_bio = get_sam_bio()

    # aggregated table of alterations fed to the heatmap
    df_sam_alt = get_sam_alt(args.sam_alt)

    # MC3 public
    df_sam_mc3_pub = get_sam_mc3_pub(args.mc3_pub)
    df_sam_mc3_con = get_sam_mc3_con(args.mc3_con)

    # join all samples tables
    df_res = outer_join_sam_tables(df_sam_awg, df_sam_bio, df_sam_alt, df_sam_mc3_pub, df_sam_mc3_con)

    # mutations for TCGA from meta-prism files
    df_mut_bio_fil, df_mut_bio_pas = get_mut_bio()

    # add info mutations from Bio PASS
    cols_show = ["HGVSc", "HGVSp_Short", "t_depth", "t_vaf", "t_ref_count", "t_alt_count", "n_depth"]
    df_res = add_gene_info_mut(df_res, df_mut=df_mut_bio_pas, col_id=col_id, gene="KRAS", suffix="Bio_PASS",
                               cols_show=cols_show)


    # mutations from MC3 controlled
    df_mut_mc3_con = get_mut_mc3_con(args.mc3_con)
    df_mut_mc3_con["FILTER"] = df_mut_mc3_con["FILTER"].apply(lambda x: x.replace("|",","))
    df_mut_mc3_con["CENTERS"] = df_mut_mc3_con["CENTERS"].apply(lambda x: x.replace("|",","))
    df_mut_mc3_con["DBVS"] = df_mut_mc3_con["DBVS"].apply(lambda x: x.replace("|",","))

    df_mut_bio_fil = df_mut_bio_fil.rename(columns={"CHROM": "Chromosome",
                                                    "POS": "Start_Position",
                                                    "Tumor_Sample": "Tumor_Sample_Barcode",
                                                    "Normal_Sample": "Matched_Norm_Sample_Barcode",
                                                    "FILTER": "FILTER_METAPRISM"})

    cols_on = ["Chromosome", "Start_Position", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"]
    df_mut_mc3_con = df_mut_mc3_con.merge(df_mut_bio_fil[cols_on+["FILTER_METAPRISM"]], how="left", on=cols_on)

    cols_show = ["HGVSc", "HGVSp_Short", "t_depth", "t_vaf", "t_ref_count", "t_alt_count", "n_depth",
                 "CENTERS", "DBVS", "FILTER", "FILTER_METAPRISM"]
    df_res = add_gene_info_mut(df_res, df_mut=df_mut_mc3_con, col_id=col_id, gene="KRAS", suffix="MC3_Con",
                               cols_show=cols_show)

    # save
    df_res.to_excel(args.output, index=False)
    print("-file saved at %s" %  args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform a detail analysis of KRAS mutations in PAAD of TCGA.")
    parser.add_argument('--sam_awg', type=str, help='Path to table of TCGA samples analyzed in the TCGA PAAD AWG.',
                        default="resources/kras_checks/NIHMS898701-supplement-2.xlsx")
    parser.add_argument('--sam_alt', type=str, help='Path to table of TCGA samples analyzed in the heatmap.',
                        default="../../../results/somatic_mutations/selection/selection_samples_tcga.tsv")
    parser.add_argument('--mc3_pub', type=str, help='Path to MC3 public PAAD table.',
                        default="../../../results/somatic_mutations/kras_checks/mc3.v0.2.8.PUBLIC_PAAD.maf.gz")
    parser.add_argument('--mc3_con', type=str, help='Path to MC3 controlled PAAD table.',
                        default="../../../results/somatic_mutations/kras_checks/mc3.v0.2.8.CONTROLLED_PAAD.maf.gz")
    parser.add_argument('--output', type=str, help='Path to output plot.',
        default="../../../results/somatic_mutations/kras_checks/results_comparison.xlsx")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
