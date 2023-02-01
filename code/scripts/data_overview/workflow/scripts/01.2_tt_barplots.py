# -*- coding: utf-8 -*-
"""
@created: Sep 27 2021
@modified: Dec 30 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Barplots showing the distribution of available data per tumor type and the distribution of histological subtypes in
rare tumors.
"""

import argparse
import pandas as pd
import sys

# local functions
sys.path.append("workflow/functions")
from util import select_tumor_types, get_table_counts
from load import load_cln_cohorts, _shorten

# pyprism, pyprismtools
from pyprism.data import load_colors
from pyprismtools import draw_barplot_counts, draw_barplot_counts_group_stack

# plotly
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pio.templates.default = "plotly_white"

def get_tumor_types_order(df):
    mask_dna = df["Sample_Type"].apply(lambda x: "DNA" in x if type(x)==str else False)
    tt_keep = df["Project_TCGA_More"].unique().tolist()

    df_dna = df.loc[mask_dna].value_counts("Project_TCGA_More")
    tt_keep_in = set(tt_keep).intersection(set(df_dna.index))
    tt_keep_other = set(tt_keep).difference(tt_keep_in)
    df_dna = df_dna.loc[tt_keep_in]
    tt_keep_order = df_dna.sort_values(ascending=False).index.tolist()

    return tt_keep_order + [x for x in tt_keep if x in tt_keep_other]

def prepare_counts_sample_type(df, tt_all=None, barmode="stack"):
    x_group = "Project_TCGA_More"
    x_stack = "Sample_Type_Simple"

    def _sample_type_simple(x):
        if type(x)==float:
            return x
        else:
            types = []
            if "DNA" in x or ("DNA" in x and "RNA" in x):
                types.append("WES or WES & RNAseq")
            elif "RNA" in x:
                types.append("RNAseq only")
            return " & ".join(types)

    df["Sample_Type_Simple"] = df["Sample_Type"].apply(_sample_type_simple)
    df_count = df.groupby([x_group, x_stack]).size().to_frame("Count").reset_index()

    stacks = ["WES or WES & RNAseq", "RNAseq only"]

    dfs_counts = {}
    for stack in stacks:
        df_stack = df_count.loc[df_count[x_stack]==stack, [x_group, "Count"]]
        df_stack = df_stack.rename(columns={x_group: "Label"})
        if tt_all is not None:
            for tt in tt_all:
                if tt not in df_stack["Label"].values:
                    row = {"Label": tt, "Count": 0}
                    df_stack = df_stack.append(row, ignore_index=True)
        df_stack = df_stack.sort_values(by="Label").reset_index(drop=True)
        dfs_counts[stack] = df_stack

    # text
    for stack in stacks:
        if stack == "WES or WES & RNAseq":
            dfs_counts[stack]["Text"] = dfs_counts[stack]["Count"].astype(str)
        elif stack == "RNAseq":
            df = dfs_counts[stack].copy()
            df = df.merge(dfs_counts["WES or WES & RNAseq"], how="left", on="Label")
            df = df.fillna(0)
            if barmode=="stack":
                df["Text"] = (df["Count_x"] + df["Count_y"]).astype(int).astype(str)
            elif barmode=="group":
                df["Text"] = (df["Count_x"]).astype(int).astype(str)
            dfs_counts[stack] = dfs_counts[stack].merge(df[["Label", "Text"]], how="left", on="Label")
    return dfs_counts


def prepare_counts_not_tcga(df, x="Histological_Type_Simple", y="Primary_Site", c="Project_TCGA"):
    yvalues = df[y].dropna().unique()
    dfs_not_tcga = {yval: df.loc[(df[c]=="Not_TCGA") & (df[y]==yval)] for yval in yvalues}
    dfs_not_tcga = {k: get_table_counts(df, x, add_pct=True) for k, df in dfs_not_tcga.items()}
    return dfs_not_tcga


def draw_barplot_counts_sample_type(dfs_counts_st, cohort, cohorts2colors, tt_keep_order, yaxis_range, barmode):
    names2colors = {"WES or WES & RNAseq": cohorts2colors["%s_dark" % cohort],
                    "RNAseq only": cohorts2colors["%s_clear" % cohort]}

    fig = draw_barplot_counts(dfs_counts_st, x="Label", y="Count", names2colors=names2colors, sort=False,
                              yaxis_range=yaxis_range, title=None, barmode=barmode, name_upper=False)
    fig.update_xaxes(tickangle=45, categoryorder='array', categoryarray=tt_keep_order)
    fig.update_layout(showlegend=True)

    return fig


def draw_group_stack_bar_tumor_types_sample_types(dfs_cln, cohorts, cohorts2colors, tt_keep_order):
    dfs_tt_st = {}
    for cohort in cohorts:
        df = dfs_cln[cohort].loc[dfs_cln[cohort]["Project_TCGA_More"].isin(tt_keep_order)].copy()
        dt = prepare_counts_sample_type(df, tt_all=tt_keep_order)

        for st in dt:
            if st == "RNAseq only":
                df_count = df.groupby("Project_TCGA_More").size().to_frame("Text").astype(str).reset_index()
                df_count["Text"] = df_count["Text"].replace("0", "")
                df_count = df_count.rename(columns={"Project_TCGA_More": "Label"})
                df_st = dt[st]
                df_st = df_st.merge(df_count, how="left", on="Label")
                dt[st] = df_st

        dfs_tt_st[cohort] = dt

    names_stacks2colors = {(c, "WES or WES & RNAseq"): cohorts2colors["%s_dark" % c] for c in cohorts}
    names_stacks2colors = {**names_stacks2colors, **{(c, "RNAseq only"): cohorts2colors["%s_clear" % c] for c in cohorts}}

    fig = draw_barplot_counts_group_stack(dfs_tt_st, names_stacks2colors, yaxis_range=[0,3.5],
                                          one_legend_per_stack=False, showlegend=True,
                                          textposition_0="outside",
                                          textfont_0=dict(color="black", family="Helvetica", size=12),
                                          textposition_1="inside",
                                          textfont_1=dict(color="white", family="Helvetica", size=12))
    fig.update_xaxes(tickangle=45, categoryorder='array', categoryarray=tt_keep_order)
    fig.update_yaxes(type="log")
    return fig


def draw_barplot_counts_not_tcga(df_cln, cohort, PSS2colors):
    x = "Histological_Type_Simple"
    y = "Primary_Site"
    dfs_counts_not_tcga = prepare_counts_not_tcga(df_cln, x, y)
    y_max = max([df["Count"].sum() for df in dfs_counts_not_tcga.values()])
    df_cln_not_tcga = df_cln.loc[df_cln["Project_TCGA"]=="Not_TCGA"]

    fig = draw_barplot_counts(dfs_counts_not_tcga, y=x, x="Count", names2colors=PSS2colors, sort=True,
                              title="Rare (Not TCGA) histological types - %s" % cohort.upper(),
                              barmode="stack", textangle=0, orientation="h")
    fig.update_yaxes(categoryorder="array", categoryarray=sorted(df_cln_not_tcga[x].unique())[::-1])
    fig.update_layout(legend=dict(font=dict(size=7), x=0.1, y=1, bgcolor='rgba(0,0,0,0)'), legend_orientation="v")
    return fig


def main(args):
    # data
    dfs_cln = load_cln_cohorts(args.cohorts)

    # replace all not Not_TCGA by rare not tcga
    type_rare = "Rare subtypes"
    for cohort, df_cln in dfs_cln.items():
        mask_rare = df_cln["Project_TCGA_More"].apply(lambda x: "Not_TCGA" in x)
        df_cln.loc[mask_rare, "Project_TCGA_More"] = type_rare
        df_cln["Project_TCGA_More"] = df_cln["Project_TCGA_More"].replace({"Unknown_Primary": "Unknown primary"})
        dfs_cln[cohort] = df_cln

    # colors
    cohorts2colors = load_colors(sheet="Global")
    PTM2colors = load_colors(sheet="Project_TCGA_More")
    PS2colors = load_colors(sheet="Primary_Site")
    PSS2colors = {_shorten(k, size=23): v for k,v in PS2colors.items()}

    # plots rare (Not_TCGA) histological_types - primary_site per cohort
    for cohort, output in zip(args.cohorts, args.outputs_rare):
        fig = draw_barplot_counts_not_tcga(dfs_cln[cohort], cohort, PS2colors)
        fig.update_layout(uniformtext_minsize=10, uniformtext_mode='show', font=dict(family="Helvetica", color="black"))
        fig.update_yaxes(tickfont=dict(size=10))
        fig.update_xaxes(tickfont=dict(size=16), showline=True, visible=True, gridwidth=2, linewidth=1,
                         linecolor="black")
        fig.update_layout(legend=dict(x=1, y=1))
        fig.write_image(output, width=0.8*args.width, height=args.height*1.75, engine="kaleido")
        print("-file saved at %s" % output)

        if cohort=="prism":
            fig.write_image(args.outputs_paper[0], width=0.8*args.width, height=args.height*1.75, engine="kaleido")
            print("-file saved at %s" % args.outputs_paper[0])

    # select tumor types
    dfs_cln["prism"] = select_tumor_types(df_cln=dfs_cln["prism"], min_count=args.min_count,
                                          tumor_types_drop=args.tt_drop)
    for cohort, df_cln in dfs_cln.items():
        if cohort!="prism":
            df_cln = select_tumor_types(df_cln=df_cln, tumor_types_keep=dfs_cln["prism"]["Project_TCGA_More"].unique())
            dfs_cln[cohort] = df_cln


    # plots tumor_types - sample_types per cohort
    for (cohort, df_cln), output, output_paper in zip(dfs_cln.items(), args.outputs_samples, args.outputs_paper[1:4]):
        if cohort=="tcga" or cohort=="met500":
            yaxis_max = round(df_cln["Project_TCGA_More"].value_counts().max()*1.2,0)
        else:
            yaxis_max = round(df_cln["Project_TCGA_More"].value_counts().max()*1.1,0)
        tt_keep_order = get_tumor_types_order(df_cln)
        type_unknown = [x for x in tt_keep_order if "Unknown" in x]
        type_last =  [type_rare] + type_unknown
        tt_keep_order = [x for x in tt_keep_order if x not in type_last] + type_last
        dfs_counts_st = prepare_counts_sample_type(df=df_cln, barmode="group")
        fig = draw_barplot_counts_sample_type(dfs_counts_st, cohort, cohorts2colors, tt_keep_order, [0,yaxis_max],
                                              barmode="group")
        fig.update_layout(uniformtext_minsize=11, uniformtext_mode='show', font=dict(family="Helvetica", color="black"))
        fig.update_yaxes(showgrid=False)
        if cohort == "met500":
            # fig.update_layout(legend=dict(xanchor="left"))
            fig.update_layout(legend=dict(xanchor="center", x=0.5, y=1))
        else:
            fig.update_layout(legend=dict(xanchor="center", x=0.5, y=1))
        fig.update_layout(xaxis=dict(tickfont=dict(size=11)), yaxis=dict(tickfont=dict(size=11)))
        if cohort=="prism":
            plot_width = 0.8*args.width
            plot_height= 0.6*args.height
        elif cohort=="met500":
            plot_width = 0.475*args.width
            plot_height= 0.8*args.height
        elif cohort=="tcga":
            plot_width = 0.5*args.width
            plot_height= 0.8*args.height
        fig.write_image(output, width=plot_width, height=plot_height, engine="kaleido")
        fig.write_image(output_paper, width=plot_width, height=plot_height, engine="kaleido")
        print("-file saved at %s" % output)
        print("-file saved at %s" % output_paper)

    # plots tumor_types - sample_types all cohorts
    tt_keep_order = get_tumor_types_order(dfs_cln["prism"])
    fig = draw_group_stack_bar_tumor_types_sample_types(dfs_cln, args.cohorts, cohorts2colors, tt_keep_order)
    fig.update_layout(uniformtext_minsize=12, uniformtext_mode='show', font=dict(family="Helvetica", color="black"))
    fig.update_layout(xaxis=dict(tickfont=dict(size=16)), yaxis=dict(tickfont=dict(size=16)))
    fig.write_image(args.output_all, width=args.width, height=args.height, engine="kaleido")
    print("-file saved at %s" % args.output_all)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohorts', type=str, nargs="+", help='Names of the cohorts.',
                        default=["prism", "met500", "tcga"])
    parser.add_argument('--min_count', type=int, default=10,
                        help="Minimum number of tumors for a tumor type to be in restricted plots.")
    parser.add_argument('--tt_drop', type=str, nargs="*", default=["N/A", "MISC - Not_TCGA"],
                        help="Tumor types excluded from plots.")
    parser.add_argument('--width', type=int, help='Width of the plot in pixels', default=1000)
    parser.add_argument('--height', type=int, help='Height of the plot in pixels', default=500)
    parser.add_argument('--outputs_samples', type=str, nargs="+", help='Paths to output plots.',
        default=["../../../results/data_overview/tumor_types/bar_chart_tumor_types_sample_types_prism.svg",
        "../../../results/data_overview/tumor_types/bar_chart_tumor_types_sample_types_met500.svg",
        "../../../results/data_overview/tumor_types/bar_chart_tumor_types_sample_types_tcga.svg"])
    parser.add_argument('--output_all', type=str, help='Paths to output plot combining all cohorts.',
        default="../../../results/data_overview/tumor_types/bar_chart_tumor_types_sample_types_all_cohorts.svg")
    parser.add_argument('--outputs_rare', type=str, nargs="+", help='Paths to output plot combining all cohorts.',
        default=["../../../results/data_overview/tumor_types/bar_chart_rare_histological_types_prism.svg",
        "../../../results/data_overview/tumor_types/bar_chart_rare_histological_types_met500.svg"])
    parser.add_argument('--outputs_paper', type=str, nargs="+", help='Paths to output plot combining all cohorts.',
        default=["../../../results/figures_paper/FS1.pdf",
                 "../../../results/figures_paper/F1b.svg",
                 "../../../results/figures_paper/FS3d.svg",
                 "../../../results/figures_paper/FS3b.svg"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
