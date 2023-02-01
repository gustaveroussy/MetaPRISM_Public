# -*- coding: utf-8 -*-
"""
@created: Jun 08 2022
@modified: Jan 02 2023
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Donut and bar charts detailing frequency of pathogenic germline mutations.
"""

import argparse
import os
import pandas as pd
import sys

# local functions
sys.path.append("../../data_overview/workflow/functions")
from load import load_cln_cohorts
from util import draw_plot_donut_count, plot_donut_count, annotation_centered_name

# pyprism
from pyprism.data import load_colors, load_table
from pyprismtools import draw_barplot_counts

# plotly
import plotly.io as pio
from plotly.subplots import make_subplots
pio.templates.default = "plotly_white"

# functions ============================================================================================================

def main(args):
    # load mutations
    dfs_mut = {c: pd.read_table(mut) for c, mut in zip(args.cohorts, args.mutations)}
    dfs_sam = {c: pd.read_table(sam) for c, sam in zip(args.cohorts, args.samples)}

    for c in args.cohorts:
        df_mut = dfs_mut[c]
        df_sam = dfs_sam[c]
        df_mut = df_mut.loc[df_mut["Sample_Id"].isin(df_sam["Sample_Id"])]
        df_mut = df_mut.merge(df_sam, how="left", on="Sample_Id")
        dfs_mut[c] = df_mut

    # load pathways
    df_mut = dfs_mut["prism"]
    df_pth = pd.read_table(args.pathways)
    df_pth = df_pth.rename(columns={"gene": "Hugo_Symbol"})
    df_mut = df_mut.merge(df_pth, how="left", on="Hugo_Symbol")

    # get count table for drawing pie charts
    dfs_dat = {}

    cohorts = [c.upper() for c in args.cohorts]
    cohorts.append("ExAC")
    fractions = [dfs_mut[c]["Subject_Id"].nunique()/dfs_sam[c]["Subject_Id"].nunique() for c in args.cohorts]
    fractions.append(0.07681093)
    df_a = pd.DataFrame({"Cohort": cohorts, "Fraction_of_pathogenic_variants": fractions})
    dfs_dat["Barplot_all"] = df_a

    df_b = df_mut["type1"].value_counts().to_frame("Count")
    df_b.index.name = "Pathway"
    df_b = df_b.reset_index(drop=False)
    dfs_dat["Piechart_prism"] = df_b

    df_mut_dna_repair = df_mut.loc[df_mut["type1"].isin(["DNA_repair"])]
    df_c = df_mut_dna_repair["type2"].value_counts().to_frame("Count")
    df_c.index.name = "Pathway"
    df_c = df_c.reset_index(drop=False)
    dfs_dat["Piechart_DNA_repair"] = df_c

    df_d = df_mut_dna_repair.groupby(["Tumor_Type", "type2"]).size().to_frame("Count")
    df_d = df_d.reset_index(drop=False).rename(columns={"type2": "Pathway"})
    dfs_dat["Piechart_DNA_repair_by_cohort"] = df_d

    # colors
    PTM2colors = load_colors(sheet="Project_TCGA_More")
    cohort2colors = load_colors(sheet="Global")
    cohort2colors = {k.upper(): v for k,v in cohort2colors.items()}
    cohort2colors["ExAC"] = "gold"
    pathways2colors = {"Cell cycle": "#b8b8ff", "Collagen": "#ffd8be", "DNA repair": "#d4e09b",
                       "Genome maintainence": "#ffa69e", "Metabolism": "#70d6ff",  "Signaling": "#ffd670"}
    dnarepair2colors = {"BER": "#d4e09b" , "FA": "#f6f4d2", "GM": "#cbdfbd", "HR": "#f19c79",
                        "MMR": "#a44a3f", "NER": "#ffddd2"}

    # figures
    figs = []

    # first plot: barplot
    name = "Barplot_all"
    dfs = {}
    for cohort, fraction in zip(dfs_dat[name]["Cohort"], dfs_dat[name]["Fraction_of_pathogenic_variants"]):
        dfs[cohort] = pd.DataFrame({"Cohort": [cohort], "Count": [fraction], "Text": ["%.3f" % fraction]})
    fig = draw_barplot_counts(dfs=dfs, x="Cohort", y="Count", names2colors=cohort2colors, textangle=0,
                              marker_line=dict(width=0.35, color='black'),
                              textfont=dict(size=12, family="Helvetica", color="black"))
    fig.update_layout(showlegend=False,
                      yaxis=dict(range=[0,0.175], tickmode='array', tickvals=[0,0.05,0.1], ticktext=["0","0.05","0.1"]))
    fig.update_yaxes(tickfont=dict(size=12, family="Helvetica", color="black"))
    fig.update_xaxes(tickangle=-90, tickfont=dict(size=12, family="Helvetica", color="black"))
    figs.append(fig)

    # second plot: donut plot
    name = "Piechart_prism"
    df = dfs_dat[name].rename(columns={"Pathway": "Label"})
    df["Label"] = df["Label"].apply(lambda x: x.replace("_", " "))
    df["Color"] = df["Label"].map(pathways2colors)
    fig = draw_plot_donut_count(cohort="PRISM", df_counts=df, size=9, size_annotation=16, textposition="outside",
                               rotation=-45, hole=0.5)
    figs.append(fig)

    # third plot: donut plot, small
    name = "Piechart_DNA_repair"
    df = dfs_dat[name].rename(columns={"Pathway": "Label"})
    df["Label"] = df["Label"].apply(lambda x: x.replace("_", " "))
    df["Color"] = df["Label"].map(dnarepair2colors)
    fig = draw_plot_donut_count(cohort="", df_counts=df, size=12, size_annotation=18, textposition="inside")
    figs.append(fig)

    # fourth plot: donut plots per tumor type
    n_row = 2
    n_col = 3
    fig = make_subplots(n_row, n_col, specs=[[{'type':'domain'}]*n_col]*n_row, horizontal_spacing=0.1,
                       vertical_spacing=0.1)
    annotations = []
    tumor_types = ["BLCA", "BRCA", "COAD", "LUAD", "LUSC", "PAAD"]
    name = "Piechart_DNA_repair_by_cohort"
    df = dfs_dat[name].rename(columns={"Pathway": "Label"})
    df["Color"] = df["Label"].map(dnarepair2colors)
    for i, tumor_type in enumerate(tumor_types):
        df_tumor_type = df.loc[df["Tumor_Type"]==tumor_type]
        trace = plot_donut_count(df_tumor_type, name=tumor_type, size=12,textposition="inside", hole=0.55)
        i_row = i // n_col + 1
        i_col = i % n_col + 1
        fig.add_trace(trace, i_row,  i_col)
        annotation = annotation_centered_name(fig_data=fig.to_dict()["data"][-1], size=18)
        annotations.append(annotation)
    fig.update_layout(annotations=annotations, showlegend=False, margin=dict(t=10, l=0, b=10, r=0))
    figs.append(fig)

    # plot
    for fig, output, output_paper, width, height in zip(figs, args.outputs, args.outputs_paper, args.widths, args.heights):
        fig.write_image(output, width=width, height=height, engine="kaleido")
        fig.write_image(output_paper, width=width, height=height, engine="kaleido")
        print("-file saved at %s" % output)
        print("-file saved at %s" % output_paper)

# Run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plots for the analysis of germline mutations.')
    parser.add_argument('--cohorts', type=str, nargs="+", help='Cohort names.',
                        default=["prism", "met500", "tcga"])
    parser.add_argument('--samples', type=str, nargs="+", help='Path to table of sample selection.',
                        default=["../../../results/germline_mutations/selection/selection_samples_prism.tsv",
                                  "../../../results/germline_mutations/selection/selection_samples_met500.tsv",
                                  "../../../results/germline_mutations/selection/selection_samples_tcga.tsv"])
    parser.add_argument('--mutations', type=str, nargs="+", help='Path table of germline mutations.',
                        default=["../../../data/prism/wes/germline_maf/germline_calls_pathogenic.maf.gz",
                                  "../../../data/met500/wes/germline_maf/germline_calls_pathogenic.maf.gz",
                                  "../../../data/tcga/wes/germline_maf/germline_calls_pathogenic.maf.gz"])
    parser.add_argument('--pathways', type=str, help='Path to table of pathways.',
                        default="resources/germline_pathways.txt")
    parser.add_argument('--widths', type=int, nargs="+", help='Widths of the plots in pixels',
                        default=[300, 250, 125, 500])
    parser.add_argument('--heights', type=int, nargs="+", help='Heights of the plots in pixels',
                        default=[300, 250, 125, 300])
    parser.add_argument('--outputs', type=str, nargs="+", help='Paths to output plots.',
                        default=["../../../results/germline_mutations/bars_pies/bars_fraction_all.pdf",
                                 "../../../results/germline_mutations/bars_pies/donut_prism_all.pdf",
                                 "../../../results/germline_mutations/bars_pies/donut_prism_all_dna_repair.pdf",
                                 "../../../results/germline_mutations/bars_pies/donut_prism_per_tumor_type.pdf"])
    parser.add_argument('--outputs_paper', type=str, nargs="+", help='Paths to output plots.',
                        default=["../../../results/figures_paper/F3c.svg",
                                 "../../../results/figures_paper/F3d_top.svg",
                                 "../../../results/figures_paper/F3d_bot.svg",
                                 "../../../results/figures_paper/F3e.svg"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
