# -*- coding: utf-8 -*-
"""
@created: Oct 15 2021
@modified: Oct 15 2021
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
from load import load_cln_cohorts
from util import get_table_counts

# pyprism
from pyprism.data import load_colors, load_table
from pyprism.util import explode_df

# pyprismtools
from pyprismtools import get_dict_labels_to_colors
from pyprismtools import draw_barplot_counts, draw_barplot_counts_group_stack

# plotly
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pio.templates.default = "plotly_white"

def prepare_counts_stack(df, x_group="Drug_Name", x_stack="Histological_Type_Simple"):
    stacks = df[x_stack].dropna().unique()
    dfs = {stack: df.loc[df[x_stack]==stack] for stack in stacks}
    dfs = {k: get_table_counts(df, x_group, add_pct=True) for k, df in dfs.items()}
    return dfs


def simplify_histological_type(x):
    neuroendocrines = ["Large cell neuroendocrine carcinoma", "Small cell neuroendocrine carcinoma",
                       "Neuroendocrine tumor"]
    if x in neuroendocrines:
        return "Neuroendocrine tumor/carcinoma"
    else:
        return x


def add_drug_class(df_drug_table, max_level=None):
    levels = sorted([x for x in df_drug_table.columns if "Class_Lvl" in x])
    if max_level is not None:
        levels = levels[:max_level]
    df_drug_table["Class"] = df_drug_table[levels].apply(lambda x: " - ".join(x.dropna().tolist()), axis=1)
    return df_drug_table


def select_sample_ids(df_alt, df_cln):
    sample_ids_dna_t = df_cln["Sample_Id_DNA_T"].dropna().tolist()
    sample_ids_rna_t = df_cln["Sample_Id_RNA_T"].dropna().tolist()
    sample_ids = sample_ids_dna_t + sample_ids_rna_t
    return df_alt.loc[df_alt["Sample_Id"].isin(sample_ids)]


def main(args):
    # load cln table
    df_cln = load_cln_cohorts([args.cohort])
    df_cln = df_cln.loc[df_cln["Project_TCGA_More"].isin(args.tumor_type.split("_"))]
    df_cln["Histological_Type_Simple"] = df_cln["Histological_Type_Simple"].apply(simplify_histological_type)

    # remove samples that have not received any drugs
    df_cln = df_cln.loc[df_cln["Drugs_Before_Biopsy"]!="None"]

    # load alt table
    df_alt = load_table(args.alts_table)
    df_alt = select_sample_ids(df_alt, df_cln)
    df_alt = df_alt.replace("", np.nan)

    # build dataframe of interest
    cols_drug = ["Subject_Id", "Project_TCGA_More", "Biopsy_Selected", "Histological_Type_Simple", "Drugs_Before_Biopsy"]
    df_drug = df_cln[cols_drug].dropna(subset=["Drugs_Before_Biopsy"])
    df_drug = df_drug.rename(columns={"Drugs_Before_Biopsy": "Drug_Name"})

    # add resistances from molecular profiles
    col = "Res_Drug_Simple"
    df_alt_res = df_alt.groupby("Subject_Id").agg({col: lambda x: "|".join(x.dropna())}).reset_index()
    df_alt_res[col] = df_alt_res[col].apply(lambda x: "|".join(sorted(list(set(x.split("|"))))))
    df_alt_res[col] = df_alt_res[col].replace("", np.nan)
    df_drug = df_drug.merge(df_alt_res, how="left", on="Subject_Id")

    # explode
    label_no_resist = "No predicted resistance"
    df_drug = explode_df(df=df_drug, cols=["Drug_Name"], sep="|")
    df_drug["Res_Drug_Simple"] = df_drug["Res_Drug_Simple"].replace(np.nan, label_no_resist)
    df_drug = explode_df(df=df_drug, cols=["Res_Drug_Simple"], sep="|")
    df_drug = df_drug.rename(columns={"Res_Drug_Simple": "Res_Drug_Name"})

    # add drug classes
    df_drug_table = load_table(args.drug_table)
    df_drug_table = add_drug_class(df_drug_table)
    dci_2_class = {r["DCI"]: r["Class"] for _,r in df_drug_table.iterrows()}
    df_drug["Drug_Class"] = df_drug["Drug_Name"].map(dci_2_class)
    df_drug["Res_Drug_Class"] = df_drug["Res_Drug_Name"].map(dci_2_class)
    df_drug.loc[df_drug["Res_Drug_Name"]==label_no_resist, "Res_Drug_Class"] = label_no_resist

    # Colors
    PTM2colors = load_colors(sheet="Project_TCGA_More")
    HTS_unique = set(df_drug["Histological_Type_Simple"])
    HTS2colors = get_dict_labels_to_colors(labels_unique=HTS_unique, sort=True, cmap_name="tab20")
    RDN_unique = set(df_drug["Res_Drug_Name"].dropna())
    RDN2colors = get_dict_labels_to_colors(labels_unique=RDN_unique, sort=True, cmap_name="tab20")
    RDN2colors[label_no_resist] = "lightgray"
    RDC_unique = set(df_drug["Res_Drug_Class"].dropna())
    RDC2colors = get_dict_labels_to_colors(labels_unique=RDC_unique, sort=True, cmap_name="tab20")
    RDC2colors[label_no_resist] = "lightgray"

    # Stacked barplot drug names | drug classes - histological type
    x_groups = ["Drug_Name", "Drug_Class"]
    x_stack = "Histological_Type_Simple"

    for x_group, width, height, output in zip(x_groups, args.widths, args.heights, args.output_histologies):
        title = "Drug resistances and histological types - %s %s" % (args.tumor_type, args.cohort.upper())
        stacks2colors =  HTS2colors
        cols = ["Subject_Id", x_group, x_stack]
        df_drug_sub = df_drug[cols].drop_duplicates()

        yaxis_range = [0, round(df_drug_sub[x_group].value_counts().max()*1.1,0)]
        x_group_order = df_drug_sub[x_group].value_counts().index.tolist()

        dfs_counts = prepare_counts_stack(df_drug_sub, x_group=x_group, x_stack=x_stack)
        fig = draw_barplot_counts(dfs_counts, x=x_group, y="Count", names2colors=stacks2colors, sort=True,
                              yaxis_range=yaxis_range, title=title, barmode="stack")
        fig.update_xaxes(tickangle=45, tickfont=dict(size=14), categoryorder='array', categoryarray=x_group_order)
        fig.update_layout(showlegend=True, legend=dict(font=dict(size=14)))
        fig.write_image(output, width=width, height=height, engine="kaleido")
        print("-file saved at %s" % output)


    # Stacked barplot drug names | drug classes - res_drug names | res_drug classes
    x_groups = ["Drug_Name", "Drug_Class"]
    x_stacks = ["Res_Drug_Name", "Res_Drug_Class"]
    stacks2colors = [RDN2colors, RDC2colors]

    for x_group, x_stack, names2colors, width, height, output in zip(x_groups, x_stacks, stacks2colors,
                                                       args.widths, args.heights, args.output_alterations):
        title = "Drug resistances and predicted resistances - %s %s" % (args.tumor_type, args.cohort.upper())
        stacks2colors =  HTS2colors
        cols = ["Subject_Id", x_group, x_stack]
        df_drug_sub = df_drug[cols].drop_duplicates()

        yaxis_range = [0, round(df_drug_sub[x_group].value_counts().max()*1.1,0)]
        x_group_order = df_drug_sub[x_group].value_counts().index.tolist()

        dfs_counts = prepare_counts_stack(df_drug_sub, x_group=x_group, x_stack=x_stack)
        fig = draw_barplot_counts(dfs_counts, x=x_group, y="Count", names2colors=names2colors, sort=True,
                              yaxis_range=yaxis_range, title=title, barmode="stack")
        fig.update_xaxes(tickangle=45, tickfont=dict(size=14), categoryorder='array', categoryarray=x_group_order)
        fig.update_layout(showlegend=True, legend=dict(font=dict(size=14)))
        fig.write_image(output, width=width, height=height, engine="kaleido")
        print("-file saved at %s" % output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Draw barplots showing drug distribution and associations with predicted resistances")
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--tumor_type', type=str, default="LUAD", help='Selection of a subset of the cohort.')
    parser.add_argument('--drug_table', type=str, help='Path to the drug table.',
                        default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
    parser.add_argument('--alts_table', type=str, help='Path to the table of annotated genomic alterations.',
                        default="../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv")
    parser.add_argument('--widths', type=int, nargs="+", help='Width of the outputs in pixels.', default=[1400, 1000])
    parser.add_argument('--heights', type=int, nargs="+", help='Height of the outputs in pixels.', default=[500, 650])
    parser.add_argument('--output_histologies', type=str, nargs="+", help='Paths to output plot with histologies.',
                    default=["../../../results/data_overview/treatments/barplots/LUAD/trt_name_vs_histologies_prism.pdf",
                             "../../../results/data_overview/treatments/barplots/LUAD/trt_class_vs_histologies_prism.pdf"])
    parser.add_argument('--output_alterations', type=str, nargs="+", help='Paths to output plot with histologies.',
                    default=["../../../results/data_overview/treatments/barplots/LUAD/trt_name_vs_alterations_prism.pdf",
                             "../../../results/data_overview/treatments/barplots/LUAD/trt_class_vs_alterations_prism.pdf"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
