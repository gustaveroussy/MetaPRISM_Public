# -*- coding: utf-8 -*-
"""
@created: 11 Jun 2022
@modified: 11 Jun 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Barplots showing the distribution of treatments (drug names and classes) across tumor types and tumor subtypes
(histological types or molecular subtypes).
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys

# local functions
sys.path.append("workflow/functions")
from util import get_table_counts

# pyprism
from pyprism.data import load_table, load_cln, load_bio, split_targeted_therapy_targets, load_colors
from pyprism.util import explode_df

# pyprismtools
from pyprismtools import draw_barplot_counts

def get_table_of_drug_data(cohort, tumor_type, level, selection={}):
    col_tt = "Project_TCGA_More"

    # load clinical and biospecimen data
    df_cln = load_cln(cohort)
    df_bio = load_bio(cohort)

    # add columns from biospecimen data
    df_cln["Biopsy_Id"] = df_cln["Biopsy_Selected"].apply(lambda x: x.split("|")[0])
    cols_bio = ["Biopsy_Id", "Biopsy_Site", "Biopsy_Subsite", "Study"]
    df_cln = df_cln.merge(df_bio[cols_bio].drop_duplicates(), how="left", on="Biopsy_Id")

    # select samples based on values if any
    if len(selection)>0:
        for col, val, mod in selection.items():
            if mod=="regex":
                r = re.compile(r'.*(%s).*' % val)
                mask = df_cln[col].apply(lambda x: bool(r.match(x)))
            elif mode=="exact":
                if type(val)!=list:
                    mask = df_cln[col]==val
                else:
                    mask = df_cln[col].isin(val)
            df_cln = df_cln[mask].copy()


    # select tumor type
    if args.tumor_type != "All":
        df_cln = df_cln.loc[df_cln[col_tt].isin(args.tumor_type.split("_"))]

    # choose column of drugs
    col_drug = "%s_Before_Biopsy" % level
    if col_drug not in df_cln.columns:
        raise ValueError("-ERROR: the level '%s' is not supported" % args.level)

    # remove patients with no data about drugs received
    mask_no_data = df_cln[col_drug].isnull()
    df_cln = df_cln.loc[~mask_no_data].copy()

    # if Classes, split by target. Drugs having multiple targets will be counted multiple times
    if level == "Classes":
        df_cln[col_drug] = df_cln[col_drug].apply(split_targeted_therapy_targets)

    # build dataframe of interest
    cols_keep = ["Subject_Id", col_tt, "Histological_Type", "Primary_Site", "Survival_Status", "Survival_Time",
                 col_drug] + cols_bio
    df_drug = df_cln[cols_keep].copy()

    # split aggregated drugs
    df_drug = explode_df(df=df_drug, cols=[col_drug], sep="|")

    return df_drug


def draw_barplot_drug_counts(df_drug, col_drug, col_tt):
    # draw simple barplot of drug counts
    x_group = col_drug
    x_stack = "Stack"

    if args.tumor_type=="All":
        df_drug[x_stack] = "All %s" % args.cohort.upper()
    else:
        df_drug[x_stack] = df_drug[col_tt]

    # get colors for stacks
    stacks2colors = load_colors(col_tt)

    # remove patients that have not received any drug
    df_drug_no_none = df_drug.loc[df_drug[col_drug]!="None"].copy()

    # get tables of percent and count for each stack
    stacks = df_drug_no_none[x_stack].dropna().unique()
    dfs_bpt = {stack: df_drug_no_none.loc[df_drug[x_stack]==stack] for stack in stacks}
    dfs_bpt = {k: get_table_counts(df, x_group, add_pct=False) for k, df in dfs_bpt.items()}
    for stack in stacks:
        df_bpt = dfs_bpt[stack]
        df_bpt["Percent"] = df_bpt["Count"]/df_drug["Subject_Id"].nunique()
        dfs_bpt[stack] = df_bpt

    # get range for y-axis
    yaxis_range = [0, round(max([df_bpt["Percent"].max() for df_bpt in dfs_bpt.values()])*1.1, 2)]
    title = "Recurrence of drugs in %s %s" % (args.tumor_type, args.cohort.upper())

    fig = draw_barplot_counts(dfs=dfs_bpt, x=x_group, y="Percent", names2colors=stacks2colors, sort=True,
                              yaxis_range=yaxis_range, title=title, barmode="stack")
    fig = fig.update_xaxes(tickangle=45, tickfont=dict(size=8))
    fig = fig.update_yaxes(title_font=dict(size=12), title_text="Percentage of patients")

    return fig


def main(args):
    col_met = "Metastatic_Sites" % args.level
    col_tt = "Project_TCGA_More"

    # load data
    df_cln = load_cln(args.cohort)
    df_bio = load_bio(cohort)

    # add columns from biospecimen data
    df_cln["Biopsy_Id"] = df_cln["Biopsy_Selected"].apply(lambda x: x.split("|")[0])
    cols_bio = ["Biopsy_Id", "Biopsy_Site", "Biopsy_Subsite", "Study"]
    df_cln = df_cln.merge(df_bio[cols_bio].drop_duplicates(), how="left", on="Biopsy_Id")

    # select samples based on values if any
    if len(selection)>0:
        for col, val, mod in selection.items():
            if mod=="regex":
                r = re.compile(r'.*(%s).*' % val)
                mask = df_cln[col].apply(lambda x: bool(r.match(x)))
            elif mode=="exact":
                if type(val)!=list:
                    mask = df_cln[col]==val
                else:
                    mask = df_cln[col].isin(val)
            df_cln = df_cln[mask].copy()


    # select tumor type
    if args.tumor_type != "All":
        df_cln = df_cln.loc[df_cln[col_tt].isin(args.tumor_type.split("_"))]

    # choose column of drugs
    col_drug = "%s_Before_Biopsy" % level
    if col_drug not in df_cln.columns:
        raise ValueError("-ERROR: the level '%s' is not supported" % args.level)

    # remove patients with no data about drugs received
    mask_no_data = df_cln[col_drug].isnull()
    df_cln = df_cln.loc[~mask_no_data].copy()

    # if Classes, split by target. Drugs having multiple targets will be counted multiple times
    if level == "Classes":
        df_cln[col_drug] = df_cln[col_drug].apply(split_targeted_therapy_targets)

    # build dataframe of interest
    cols_keep = ["Subject_Id", col_tt, "Histological_Type", "Primary_Site", "Survival_Status", "Survival_Time",
                 col_drug] + cols_bio
    df_drug = df_cln[cols_keep].copy()

    # split aggregated drugs
    df_drug = explode_df(df=df_drug, cols=[col_drug], sep="|")

    return df_drug


    # get data
    df_drug = get_table_of_drug_data(cohort=args.cohort, tumor_type=args.tumor_type, level=args.level)

    # add classes if applicable
    if args.level == "Drugs":
        cols_lvl = [x for x in df_drug_table if "_Lvl_" in x]
        df_drug_table_unique = df_drug_table[["DCI"]+cols_lvl].drop_duplicates().rename(columns={"DCI": col_drug})
        df_drug = df_drug.merge(df_drug_table_unique, how="left", on=col_drug)

    # # select drugs
    # drugs_keep = df_drug_table.loc[df_drug_table["Class_Lvl_1"]=="Chemo_Non_Alkylating"]["DCI"].unique()
    # df_drug = df_drug.loc[df_drug[col_drug].isin(drugs_keep)]

    # draw barplot showing the percentage and count of each drug
    fig = draw_barplot_drug_counts(df_drug=df_drug, col_drug=col_drug, col_tt=col_tt)
    width, height = 200 + 20*df_drug[col_drug].nunique(), 650
    fig.write_image(args.output_barplot, width=width, height=height, engine="kaleido")
    print("-file saved at %s" % args.output_barplot)

    # draw double heatmap showing cooccurrence of drugs
    fig = draw_double_heatmap_drugs(df_drug, col_drug)
    plt.subplots_adjust(bottom=0.25, left=0.25)
    fig.savefig(args.output_dbl_htp, dpi=300)
    print("-file saved at %s" % args.output_dbl_htp)

    # draw comut (oncoplot-like) plot
    plt.rcParams["legend.title_fontsize"] = 16
    plt.rcParams["legend.fontsize"] = 10
    plt.rcParams["xtick.labelsize"] = 0
    plt.rcParams["ytick.labelsize"] = 7
    plt.rcParams["axes.titlesize"] = 12

    x_padding = 0.04 # the x distance between patches in comut
    y_padding = 0.04 # the y distance between patches in comut
    tri_padding = 0.03 # the distance between triangles in comut
    figsize = (20,8)
    dpi = 300

    cols_comut = ["sample", "category", "value"]
    samples = df_drug["Subject_Id"].unique().tolist()
    df_drug_sub = df_drug.loc[df_drug[col_drug]!="None"].copy()
    drug_order = df_drug_sub[col_drug].value_counts().index.tolist()[::-1]

    # samples order using lexicographical ordering
    samples_no_drug = sorted(list(set(samples).difference(set(df_drug_sub["Subject_Id"].tolist()))))
    df_drug_sub["Received"] = 1
    df_drug_binary = df_drug_sub.pivot(index=col_drug, columns="Subject_Id", values="Received")
    df_drug_binary = df_drug_binary.loc[drug_order[::-1],:].fillna("0")
    df_drug_binary[df_drug_binary==1] = "1"
    df_drug_string = df_drug_binary.apply(lambda x: "".join(x.tolist()), axis=0)
    samples_drug_ordered = df_drug_string.sort_values(ascending=False).index.tolist()
    samples_ordered = samples_drug_ordered + samples_no_drug

    object_comut = comut.CoMut()
    object_comut.samples = samples_ordered

    # add mutations data
    df_comut = df_drug_sub.rename(columns={"Subject_Id": "sample", col_drug: "category"})
    if args.level=="Classes":
        name_drug = "Drug"
        df_comut["value"] = "Drug received"
        mapping = {"Drug received": "#F49F0A"}
    else:
        name_drug = "Drug category"
        df_comut["value"] = df_comut["Class_Lvl_1"]
        mapping = load_colors("Class_Lvl_1")

    object_comut.add_categorical_data(data=df_comut, name=name_drug, category_order=drug_order, mapping=mapping)
    # object_comut = object_comut.plot_comut(x_padding=x_padding, y_padding=y_padding, tri_padding=tri_padding,
    #                                        figsize=figsize)
    # lgd_comut = object_comut.add_unified_legend()

    # add tumor type data if relevant
    mapping = load_colors(col_tt)

    if not args.tumor_type in mapping.keys():
        name = "Tumor type"
        df_comut = df_drug.rename(columns={"Subject_Id": "sample", col_tt: "value"})
        df_comut["category"] = name
        df_comut = df_comut[cols_comut].drop_duplicates()

        object_comut.add_categorical_data(df_comut, name=name, mapping=mapping)
        object_comut = object_comut.plot_comut(x_padding=x_padding, y_padding=y_padding, tri_padding=tri_padding,
                                               hspace=0.03, figsize=figsize)
        lgd_comut = object_comut.add_unified_legend()

    # add Survival data
    name = "Survival status"
    df_comut = df_drug.rename(columns={"Subject_Id": "sample", "Survival_Status": "value"})
    df_comut["category"] = name
    df_comut = df_comut[cols_comut].drop_duplicates()
    mapping = load_colors("Survival_Status")

    object_comut.add_categorical_data(df_comut, name=name, mapping=mapping)
    # object_comut = object_comut.plot_comut(x_padding=x_padding, y_padding=y_padding, tri_padding=tri_padding,
    #                                        hspace=0.03, figsize=figsize)
    # lgd_comut = object_comut.add_unified_legend()

    # add side barplot 
    df_bar = df_drug_sub[col_drug].value_counts().reset_index()
    df_bar.columns = ["category", "count"]
    mapping = {"count": "darkgrey"}
    name = "Count"

    plt.rcParams["xtick.labelsize"] = 12
    object_comut.add_side_bar_data(df_bar, paired_name=name_drug, name=name, position='left', mapping=mapping,
                                   xlabel='# patients')
    object_comut = object_comut.plot_comut(x_padding=x_padding, y_padding=y_padding, tri_padding=tri_padding,
                                           hspace=0.03, figsize=figsize)
    lgd_comut = object_comut.add_unified_legend(axis_name=name_drug)

    object_comut.figure.savefig(args.output_comut, dpi=dpi, bbox_inches='tight')


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Draw barplots showing drug distribution.")
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--tumor_type', type=str, default="All", help='Selection of a subset of the cohort.')
    parser.add_argument('--output', type=str, help='Path to output barplot.',
                        default="../../../results/data_overview/metastatic_sites/barplot_All.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
