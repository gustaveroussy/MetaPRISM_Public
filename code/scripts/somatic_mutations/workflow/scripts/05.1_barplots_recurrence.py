# -*- coding: utf-8 -*-
"""
@created: Nov 8 2021
@modified: Feb 03 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Barplots showing the recurrence of alterations.
"""

import argparse
import numpy as np
import pandas as pd
import re
import sys
sys.path.append("../../data_overview/workflow/functions")

from util import get_table_counts

# pyprism
from pyprism.data import load_colors, load_table, preprocess_wes_mut
from pyprism.util import explode_df

# pyprismtools
from pyprismtools import get_dict_labels_to_colors
from pyprismtools import draw_barplot_counts, draw_barplot_counts_group_stack

# plotly
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pio.templates.default = "plotly_white"

# functions ============================================================================================================

def correct_dtype_from_num_to_str(x):
    try:
        y = "%d" % int(float(x))
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = x
        except:
            y = x

    return y


def add_row_identifier(df_mut):
    cols_required = ["Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                     "Chromosome", "Start_Position", "End_Position",
                     "Reference_Allele", "Tumor_Seq_Allele2"]
    df_muc = df_mut.copy()
    for col in cols_required:
        df_muc[col] = df_muc[col].apply(correct_dtype_from_num_to_str)
    df_mut["Row_Identifier"] = df_muc[cols_required].fillna("NA").astype(str).apply("-".join, axis=1)
    assert df_mut["Row_Identifier"].nunique()==df_mut.shape[0]
    return df_mut


def add_mut_identifier(df_mut):
    cols_required = ["Chromosome", "Start_Position", "End_Position",
                     "Reference_Allele", "Tumor_Seq_Allele2"]
    df_muc = df_mut.copy()
    for col in cols_required:
        df_muc[col] = df_muc[col].apply(correct_dtype_from_num_to_str)
    df_mut["Mutation_Identifier"] = df_muc[cols_required].fillna("NA").astype(str).apply("-".join, axis=1)
    return df_mut


def add_mut_identifier_label(df_mut):
    col_y = "Mutation_Identifier"
    col_x = "Mutation_Identifier_Label"
    mask_null = df_mut["HGVSp_Short"].isnull()

    cols_pc = ["Hugo_Symbol","HGVSp_Short"]
    df_mut.loc[~mask_null, col_x] = df_mut.loc[~mask_null, cols_pc[0]] + " - " + \
                df_mut.loc[~mask_null, cols_pc[1]].str.replace("%3D", "=")

    cols_pos = ["Hugo_Symbol","Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
    df_muc = df_mut.copy()
    for col in cols_pos:
        df_muc[col] = df_muc[col].apply(correct_dtype_from_num_to_str)

    df_muc.loc[mask_null, cols_pos] = df_muc.loc[mask_null, cols_pos].fillna("-")
    df_mut.loc[mask_null, col_x] = df_muc.loc[mask_null, cols_pos[0]] + " - " + \
            df_muc.loc[mask_null, cols_pos[1]] + " (" + \
            df_muc.loc[mask_null, cols_pos[2]] + ">" + \
            df_muc.loc[mask_null, cols_pos[3]] + ")"

    return df_mut


def add_variant_classification_custom(df_mut):
    old2new = {"5'UTR":"Upstream",
               "5'Flank":"Upstream",
               "3'UTR":"Downstream",
               "3'Flank":"Downstream",
               "De_novo_Start_InFrame":"Missense",
               "De_novo_Start_OutOfFrame":"Missense",
               "Frame_Shift_Del":"FSIndel",
               "Frame_Shift_Ins":"FSIndel",
               "In_Frame_Del":"IFIndel",
               "In_Frame_Ins":"IFIndel",
               "Start_Codon_Del":"StartSite",
               "Start_Codon_Ins":"StartSite",
               "Stop_Codon_Del":"Nonstop",
               "Stop_Codon_Ins":"Nonstop",
               "Start_Codon_SNP":"StartSite",
               "Translation_Start_Site":"StartSite",
               "Missense_Mutation":"Missense",
               "Nonsense_Mutation":"Nonsense",
               "Nonstop_Mutation":"Nonstop",
               "Silent":"Silent",
               "Intron":"Silent",
               "Splice_Site":"Splice",
               "Splice_Region":"Splice",
               "IGR":"IGR",
               "lincRNA": "RNA",
               "RNA":"RNA"}

    col_x = "Variant_Classification"
    col_y = "Variant_Classification_Custom"

    # convert and check that no NaN was introduced from mapping
    df_mut[col_y] = df_mut[col_x].map(old2new)
    assert df_mut[col_y].isnull().sum() == df_mut[col_x].isnull().sum()

    # add a Multihit category
    cols_id = ["Tumor_Sample_Barcode", "Hugo_Symbol"]
    df_count = df_mut.groupby(cols_id).size().to_frame("Count").reset_index()
    df_count_multi = df_count.loc[df_count["Count"] > 1]

    # remove when Hugo_Symbol="Unknown"
    df_count_multi = df_count_multi[df_count_multi["Hugo_Symbol"]!="Unknown"]
    df_count_multi["Identifier"] = df_count[cols_id].apply(" - ".join, axis=1)
    s_identifier = df_mut[cols_id].apply(" - ".join, axis=1)
    mask_multi = s_identifier.isin(df_count_multi["Identifier"])
    df_mut.loc[mask_multi, "Variant_Classification_Custom"] = "Multihit"

    return df_mut


def get_colors_x_group_x_stack(df_mut, args):
    try:
        colors = load_colors(sheet=args.color_by)
    except:
        labels_unique = df_mut[args.color_by].unique()
        if len(labels_unique)<=8:
            cmap_name="Set2"
        elif len(labels_unique)<=10:
            cmap_name="tab10"
        else:
            cmap_name="tab20"

        colors = get_dict_labels_to_colors(labels_unique=labels_unique, sort=True, cmap_name=cmap_name)

    if args.group_by=="mutations":
        x_group = "Mutation_Identifier_Label"
    elif args.group_by=="genes":
        x_group = "Hugo_Symbol"
    else:
        raise ValueError("-unsupported value %s for --group_by" % x_group)

    x_stack = args.color_by

    return colors, x_group, x_stack


def main(args):
    # load mutations and add bio, cln attributes
    df_mut_all = load_table(args.mutations_all, header_prefix="##")
    df_mut_ann = load_table(args.mutations_ann, header_prefix="##")

    # add bio, cln attributes and correct some stuff
    if args.selection_mut=="annotated":
        df_mut_ann = preprocess_wes_mut(df_mut_ann, args.cohort, cols_cln=["Project_TCGA_More"], select_pairs=True,
                                        selection_mut=args.selection_mut)
        df_mut_all = df_mut_ann
    else:
        df_mut_all = preprocess_wes_mut(df_mut_all, args.cohort, cols_cln=["Project_TCGA_More"], select_pairs=True,
                                        selection_mut=args.selection_mut)
        df_mut_ann = preprocess_wes_mut(df_mut_ann, args.cohort, cols_cln=["Project_TCGA_More"], select_pairs=True,
                                        selection_mut=args.selection_mut)

    # add identifiers
    df_mut_all = add_row_identifier(df_mut_all)
    df_mut_all = add_mut_identifier(df_mut_all)
    df_mut_ann = add_row_identifier(df_mut_ann)
    df_mut_ann = add_mut_identifier(df_mut_ann)

    # add Mutation_Identifier_Label
    df_mut_ann = add_mut_identifier_label(df_mut_ann)
    df_mut_all = add_mut_identifier_label(df_mut_all)

    # add Annotation_Status
    mask_ann = df_mut_all["Mutation_Identifier"].isin(df_mut_ann["Mutation_Identifier"])
    df_mut_all.loc[mask_ann, "Annotation_Status"] = "Annotated"
    df_mut_all.loc[~mask_ann, "Annotation_Status"] = "Not annotated"

    # add Variant_Classification_Custom
    df_mut_all = add_variant_classification_custom(df_mut_all)

    # barplot
    colors, x_group, x_stack = get_colors_x_group_x_stack(df_mut_all, args)
    title = "Recurrence of events in %s %s" % (args.tumor_type, args.cohort.upper())
    cols = ["Row_Identifier", x_group, x_stack]
    df_mut = df_mut_all[cols].drop_duplicates()

    # select recurrent events
    x_group_order = df_mut[x_group].value_counts()
    x_group_order = x_group_order[x_group_order >= args.recurrence_threshold].index.tolist()
    df_mut = df_mut.loc[df_mut[x_group].isin(x_group_order)]
    yaxis_range = [0, round(df_mut[x_group].value_counts().max()*1.1,0)]

    # counts
    stacks = df_mut[x_stack].dropna().unique()
    dfs = {stack: df_mut.loc[df_mut[x_stack]==stack] for stack in stacks}
    dfs = {k: get_table_counts(df, x_group, add_pct=True) for k, df in dfs.items()}

    fig = draw_barplot_counts(dfs, x=x_group, y="Count", names2colors=colors, sort=True,
                              yaxis_range=yaxis_range, title=title, barmode="stack")
    fig.update_xaxes(tickangle=45, tickfont=dict(size=8), categoryorder='array', categoryarray=x_group_order)
    fig.update_yaxes(title_font=dict(size=12), title_text="Number of tumor samples")

    # compute width and height
    width = 30*len(x_group_order)
    height = 650

    fig.write_image(args.output, width=width, height=height, engine="kaleido")
    print("-file saved at %s" % args.output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Draw barplots showing the recurrence of mutations.")
    parser.add_argument('--cohort', type=str, default="prism", help='Name of the cohort.')
    parser.add_argument("--mutations_all", type=str, help="Path to raw mutations table.",
                    default="../../../data/prism/wes/somatic_maf/somatic_calls.maf.gz")
    parser.add_argument("--mutations_ann", type=str, help="Path to annotated pathogenic mutations table.",
                    default="../../../data/prism/wes/somatic_maf/somatic_calls_union_ann.maf.gz")
    parser.add_argument("--selection_mut", type=str, default="non_synonymous", help="For selecting mutations.")
    parser.add_argument('--tumor_type', type=str, default="All", help='For selecting a subset of the cohort.')
    parser.add_argument('--color_by', type=str, default="Annotation_Status",
                        help='Name of the attribute used to color the bar stacks')
    parser.add_argument('--group_by', type=str, default="mutations",
                        help='Name of the attribute used to aggregate counts. You may choose "genes" or "mutations".')
    parser.add_argument('--recurrence_threshold', type=float, default=5,
                        help='Recurrence threshold for including a recurrent event in the plot.')
    parser.add_argument('--output', type=str, help='Path to output plot.',
        default="../../../results/somatic_mutations/recurrence/barplots_non_synonymous_mutations_cby_Annotation_Status_gby_genes_in_All_of_prism.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
