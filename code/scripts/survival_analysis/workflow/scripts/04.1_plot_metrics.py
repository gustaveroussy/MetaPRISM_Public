# -*- coding: utf-8 -*-
"""
@created: Dec 09 2021
@modified: May 16 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Plot c-index and brier scores estimates across different models and features selection.
"""

import argparse
import os
import numpy as np
import pandas as pd
import re
import sys
sys.path.append("workflow/functions")

from utils import load_met_cov

# pyplot
import matplotlib.cm as cm

# pyprismtools
from pyprismtools import draw_barplot_counts_group_stack

# plotly
import plotly.io as pio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
pio.templates.default = "plotly_white"

def make_name(x):
    cgby = (re.sub("coxph_", "",  x[0])).title()
    split = x[1]
    return " ".join([cgby, split])


def draw_barplots_top(df_stk, groups_ordered, dodges_ordered, x_group="Features", x_stack="Stack", x_dodge="Model",
                      **kwargs):
    x_group = "Features"
    x_stack = "Stack"
    x_dodge = "Model"

    stacks_ordered = ["No error & no warning", "Error", "Warning"]

    mask_e = df_stk["Error"].isnull()
    mask_w = df_stk["Warning"].isnull()
    stacks2colors = {"No error & no warning": "#595959", "Warning": "#E38F35", "Error": "#FF3864"}
    stacks2colors = {s: stacks2colors[s] for s in stacks_ordered}

    for stack in stacks_ordered:
        if stack == "No error & no warning":
            mask = mask_e & mask_w
        elif stack == "Error":
            mask = ~mask_e
        elif stack == "Warning":
            mask = ~mask_w
        df_stk.loc[mask, x_stack] = stack

    # dataframe for adding triplets (group, dodge, stack) with 0 count
    stacks_zer = stacks_ordered*len(groups_ordered)
    groups_zer = [g for g in groups_ordered for stack in stacks_ordered]
    df_zer = pd.DataFrame({x_stack: stacks_zer, x_group: groups_zer})

    dfs_cnt = {}
    for d in dodges_ordered:
        df_stk_d = df_stk.loc[df_stk[x_dodge]==d]
        df_zer_d = df_zer.loc[df_zer[x_group].isin(df_stk_d[x_group].unique())]
        df_cnt_d = df_stk_d.groupby([x_group, x_stack]).size().to_frame("Count").reset_index()
        df_cnt_d = df_cnt_d.merge(df_zer_d, how="outer", on=[x_stack, x_group]).fillna(0)

        dfs_cnt_d = {}
        for stk in df_stk[x_stack].unique():
            dfs_cnt_d_stk = df_cnt_d.loc[df_cnt_d[x_stack]==stk, [x_group, "Count"]]
            dfs_cnt_d_stk = dfs_cnt_d_stk.rename(columns={x_group: "Label"})
            dfs_cnt_d_stk = dfs_cnt_d_stk.sort_values(by="Label").reset_index(drop=True)
            dfs_cnt_d[stk] = dfs_cnt_d_stk

        dfs_cnt[d] = dfs_cnt_d

    names_stacks2colors = {(d, s): c for s, c in stacks2colors.items() for d in df_stk[x_dodge]}
    traces_stk = draw_barplot_counts_group_stack(dfs_cnt, names_stacks2colors, textposition_0="none",
                                                 textposition_1="none", return_traces=True, one_legend_per_stack=True,
                                                 **kwargs)
    return traces_stk


def draw_boxplot_bot(df_box, x="Features", y="C_ipcw", z="Name", z_order=None, z_colors=None, **kwargs):
    if z_order is None:
        z_order = sorted(df_box[z].unique.tolist())

    if z_colors is None:
        if len(z_order) < 10:
            cmap = "tab10"
        else:
            cmap = "tab20"
        colors_unique = cm.get_cmap(cmap).colors
        colors_unique = ["rgb(%s,%s,%s)" % color for color in colors]
        z_colors = {z_ord: color for z_ord, color in zip(z_order, colors_unique)}

    # draw series of boxplots for each z value
    data = []
    for z_ord in z_order:
        box = go.Box(
            name=z_ord,
            x=df_box.loc[df_box[z]==z_ord, x],
            y=df_box.loc[df_box[z]==z_ord, y],
            marker_color=z_colors[z_ord],
            **kwargs
        )
        data.append(box)

    return data


def main(args):
    # load pooled metrics and coefficients estimates
    df_met, df_cov = load_met_cov(dir_models=args.dir_models,
                                  dir_data=args.dir_data,
                                  names_features=args.names_features,
                                  names_models=args.names_models,
                                  names_selections=args.names_selections)

    # pool metric estimates across imputations by mean
    df_met_rep = df_met.loc[df_met["Run"]!="Full"].copy()
    del df_met

    # choose how boxplots are grouped
    if args.group_by=="Model":
        if len(args.names_selections)!=1:
            raise ValueError("-when grouping by model, you may specify only 1 selection")
        df_met_rep["Name"] = df_met_rep[["Model", "Split"]].apply(make_name, axis=1)
    elif args.group_by=="Selection":
        if len(args.names_models)!=1:
            raise ValueError("-when grouping by selection, you may specify only 1 model")
        df_met_rep["Name"] = df_met_rep[["Selection", "Split"]].apply(make_name, axis=1)

    # drop some models if asked by uesr
    mask_e = df_met_rep["Error"].isnull()
    mask_w = df_met_rep["Warning"].isnull()

    if args.drop_error=="true" and args.drop_warn=="true":
        mask_filt = mask_e & mask_w
        print("INFO: dropping %d/%d lines from models with errors or warnings." % (sum(~mask_filt), len(mask_filt)))
    elif args.drop_error=="true":
        mask_filt = mask_e
        print("INFO: dropping %d/%d lines from models with errors." % (sum(~mask_filt), len(mask_filt)))
    elif args.drop_warn=="true":
        mask_filt = mask_w
        print("INFO: dropping %d/%d lines from models with warnings." % (sum(~mask_filt), len(mask_filt)))
    else:
        mask_filt = pd.Series(True, index=df_met_rep.index)

    # draw plots comparing models using quality metrics
    colors_ordered = ["rgb(%s,%s,%s)" % color for color in cm.get_cmap("Paired").colors]
    groups_ordered = args.names_features
    dodges_ordered = args.names_models

    df_names_models = df_met_rep[["Name", "Model"]].drop_duplicates()
    df_names_models["Model"] = pd.Categorical(df_names_models["Model"], categories=dodges_ordered[::-1])
    df_names_models = df_names_models.sort_values(by=["Model", "Name"], ascending=False)

    names_ordered = df_names_models["Name"].tolist()
    names2colors = {name: color for name, color in zip(names_ordered, colors_ordered)}

    # choose time horizon for metrics
    if args.time_horizon not in df_met_rep["Time"].unique().tolist():
        preselected_times = df_met_rep["Time"].unique().tolist()
        raise ValueError("-the specified time horizon %d is not in the preselected list of time horizons:\n\t%s" %
                         (args.time_horizon, preselected_times))
    mask_time = df_met_rep["Time"]==args.time_horizon

    # choose only train splits for stack bars
    mask_train = df_met_rep["Split"]=="Train"

    # tables for box plots and stacked bars
    df_box = df_met_rep.loc[mask_filt & mask_time].copy()
    df_stk = df_met_rep.loc[mask_train & mask_time].copy()
    n_box = df_box["Features"].nunique()
    n_run = df_stk["Run"].nunique()
    width = args.width_met_one*n_box + 300

    # top stacked barplot
    traces_stk = draw_barplots_top(df_stk=df_stk, groups_ordered=groups_ordered, dodges_ordered=dodges_ordered,
                                   x_group="Features", x_stack="Stack", x_dodge="Model", xaxis="x1", yaxis="y1",
                                   legendgroup="1")

    # brier score figure
    fig = make_subplots(rows=2, cols=1, row_heights=[0.3,0.7], vertical_spacing=0.05)
    traces_met =  draw_boxplot_bot(df_box=df_box, x="Features", y="Brier", z="Name", z_order=names_ordered,
                                   z_colors=names2colors, xaxis="x2", yaxis="y2", legendgroup="2")

    for trace_stk in traces_stk:
        fig.append_trace(trace_stk, row=1, col=1)
    for trace_met in traces_met:
        fig.append_trace(trace_met, row=2, col=1)

    fig.update_layout(xaxis=dict(showticklabels=False, categoryorder='array', categoryarray=groups_ordered))
    fig.update_layout(xaxis2=dict(categoryorder='array', categoryarray=groups_ordered, title_font=dict(size=20),
                                  tickfont=dict(size=24), tickangle=45))
    fig.update_layout(yaxis_title="# of models", yaxis=dict(title_font=dict(size=24), tickfont=dict(size=18),
                                                            tickmode="array", tickvals=[0, n_run/2, n_run]))
    fig.update_layout(yaxis2_title="Brier score (ipcw)", yaxis2=dict(title_font=dict(size=24), tickfont=dict(size=18)))
    fig.update_layout(boxmode="group", legend=dict(font=dict(size=20)), legend_tracegroupgap=180)

    fig.write_image(args.output_met_bscore, width=width, height=args.height_met, engine="kaleido")
    print("-boxplot saved at %s" % args.output_met_bscore)

    # c-index figure
    fig = make_subplots(rows=2, cols=1, row_heights=[0.3,0.7], vertical_spacing=0.05)
    traces_met =  draw_boxplot_bot(df_box=df_box, x="Features", y="C_ipcw", z="Name", z_order=names_ordered,
                                   z_colors=names2colors, xaxis="x2", yaxis="y2", legendgroup="2")

    for trace_stk in traces_stk:
        fig.append_trace(trace_stk, row=1, col=1)
    for trace_met in traces_met:
        fig.append_trace(trace_met, row=2, col=1)

    fig.update_layout(xaxis=dict(showticklabels=False, categoryorder='array', categoryarray=groups_ordered))
    fig.update_layout(xaxis2=dict(categoryorder='array', categoryarray=groups_ordered, title_font=dict(size=20),
                                  tickfont=dict(size=24), tickangle=45))
    fig.update_layout(yaxis_title="# of models", yaxis=dict(title_font=dict(size=24), tickfont=dict(size=18),
                                                            tickmode="array", tickvals=[0, n_run/2, n_run]))
    fig.update_layout(yaxis2_title="C-index (ipcw)", yaxis2=dict(title_font=dict(size=24), tickfont=dict(size=18)))
    fig.update_layout(boxmode="group", legend=dict(font=dict(size=20)), legend_tracegroupgap=180)

    fig.write_image(args.output_met_cindex, width=width, height=args.height_met, engine="kaleido")
    print("-boxplot saved at %s" % args.output_met_cindex)


# Run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Draw plot showing model qualities and coefficient values.')
    parser.add_argument("--names_features", type=str, nargs="+", help="Name of the features selection.",
                        default=["cln_grim_0", "cln_biol_0", "cln_biol_1", "cln_biol_2",
                                 "cln_biol_1_dna_1", "cln_biol_1_dna_2", "cln_biol_1_dna_3", "cln_biol_1_dna_4",
                                 "cln_biol_1_dna_5", "cln_biol_1_rna_1", "cln_biol_1_rna_2",
                                 "cln_biol_1_dna_1_rna_1", "cln_biol_1_dna_5_rna_1"])
    parser.add_argument("--names_models", type=str, nargs="+", help="Name of the model.",
                        default=["coxph_standard", "coxph_ridge", "coxph_lasso"])
    parser.add_argument("--names_selections", type=str, nargs="+", help="Name of the model.",
                        default=["none"])
    parser.add_argument("--group_by", type=str,  help="Choose 'Model' or 'Selection'.",
                        default="Model")
    parser.add_argument("--dir_data", type=str, help="Path to directory of run results.",
            default="../../../results/survival_analysis/data/sub_features/BRCA_dna_and_rna")
    parser.add_argument("--dir_models", type=str, help="Path to directory of run results.",
            default="../../../results/survival_analysis/models/BRCA_dna_and_rna")
    parser.add_argument("--time_horizon", type=int, help="Time horizon for the c-index and brier score in days.",
                        default=180)
    parser.add_argument("--drop_warn", type=str, default="true",
                      help="If true, metrics from model(s) (if imputations) with warning are dropped.")
    parser.add_argument("--drop_error", type=str, default="true",
                      help="If true, metrics from model(s) (if imputations) with error are dropped.")
    parser.add_argument('--width_met_one', type=int, help='Width in pixels of 1 group of bars in boxplots.',default=140)
    parser.add_argument('--height_met', type=int, help='Height in pixels of boxplots.', default=1000)
    parser.add_argument('--output_met_bscore', type=str, help='Path to boxplots of Brier scores.',
        default="../../../results/survival_analysis/plots/BRCA_dna_and_rna/boxplot_bscore_180_none.pdf")
    parser.add_argument('--output_met_cindex', type=str, help='Path to boxplots of c-indices.',
        default="../../../results/survival_analysis/plots/BRCA_dna_and_rna/boxplot_cindex_180_none.pdf")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))

    print("\n", end="")

    main(args)
