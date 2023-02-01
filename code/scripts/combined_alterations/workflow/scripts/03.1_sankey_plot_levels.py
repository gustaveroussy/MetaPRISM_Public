# -*- coding: utf-8 -*-
"""
@created: 18 Mar 2022
@modified: 06 Jul 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Sankey plots showing the relationships between sensitivity tiers, resistance tiers and tumor type.
"""

import argparse
import pandas as pd
import sys

# pyprism
from pyprism.data import load_colors, load_table
from pyprism.util import explode_df

# holoviews with bokeh backend
import holoviews as hv
from bokeh.io import export_svgs
from bokeh.core.properties import value, field
hv.extension('bokeh')

# pyprismtools
from pyprismtools import get_dict_labels_to_colors, draw_sankey_plot

def get_level_categories_per_patient(df_alt, col_sid, sub_keep):
    dfs_alt_lvl = []
    for col_lvl in ["Sen_Level_Simple", "Res_Level_Simple"]:
        # select 1 category per patient
        df_alt_lvl = df_alt.sort_values(by=[col_sid, col_lvl])[[col_sid, col_lvl]]
        df_alt_lvl = df_alt_lvl.drop_duplicates(subset=[col_sid], keep="first")

        # add patients with no oncogenic alteration
        sub_keep_mis = list(set(sub_keep).difference(set(df_alt_lvl[col_sid].tolist())))
        df_alt_lvl_mis = pd.DataFrame({col_sid: sub_keep_mis})
        df_alt_lvl = pd.concat((df_alt_lvl, df_alt_lvl_mis))
        df_alt_lvl[col_lvl] = df_alt_lvl[col_lvl].fillna("No oncogenic/level")
        df_alt_lvl[col_lvl] = df_alt_lvl[col_lvl] + " " + col_lvl[0]

        dfs_alt_lvl.append(df_alt_lvl.set_index(col_sid))

    return pd.concat(dfs_alt_lvl, axis=1).reset_index()


def main(args):
    # load alterations table and samples selection table
    df_cln = pd.read_table(args.cln_table)
    df_alt = pd.read_table(args.alt_table)
    df_sam = pd.read_table(args.sam_table)

    mask_dna_and_rna = df_cln["Sample_Type"]=="DNA_N|DNA_T|RNA_T"
    df_counts_tt = df_cln.loc[mask_dna_and_rna]["Project_TCGA_More"].value_counts()

    # selet samples
    sam_keep = df_sam.loc[df_sam["Use_sankey"]==1]["Sample_Id"].tolist()
    sub_keep = df_cln.loc[df_cln["Sample_Id_DNA_T"].isin(sam_keep)]["Subject_Id"].tolist()
    df_alt = df_alt.loc[df_alt["Sample_Id"].isin(sam_keep)].copy()

    # categorize subjects per each of the following category
    # - has >=1 tier1 sensitivity/resistance alteration
    # - has >=1 tier2 sensitivity/resistance alteration and no tier1 sensitivity/resistance alteration
    # - has >=1 tier3 sensitivity/resistance alteration and no tier1/tier2 sensitivity/resistance alteration
    # - has no tier1/tier2/tier3 sensitivity/resistance alteration
    col_sid = "Subject_Id"
    df_alt_lvl = get_level_categories_per_patient(df_alt=df_alt, col_sid=col_sid, sub_keep=sub_keep)
    df_alt_lvl = df_alt_lvl.merge(df_cln, how="left", on="Subject_Id")

    # Select tumor types with sufficiently many tumors
    df_alt_lvl_min = df_alt_lvl.loc[df_alt_lvl["Project_TCGA_More"].isin(df_counts_tt[df_counts_tt>=args.min_count].index)]

    # remove some tumor types
    df_alt_lvl_min = df_alt_lvl_min.loc[~df_alt_lvl_min["Project_TCGA_More"].isin(["MISC - Not_TCGA"])].copy()


    # Colors
    PTM2colors = load_colors(sheet="Project_TCGA_More")
    Sen2colors = {"Tier1 S":"#00af56",
                  "Tier2 S":"#6fbfdd",
                  "Tier3 S":"#d4a8d0",
                  "No oncogenic/level S":"#ffe0c3"}

    Res2colors = {"Tier1 R":"#dd2d4a",
                  "Tier2 R":"#ffb700",
                  "Tier3 R":"#e0aaff",
                  "No oncogenic/level R":"#ffe0c3"}

    # Hook for customizing fonts
    def hook(plot, element):
        plot.handles['text_1_glyph'].text_font = value('Helvetica')
        plot.handles['text_1_glyph'].text_font_size = '12pt'
        plot.handles['text_2_glyph'].text_font = value('Helvetica')
        plot.handles['text_2_glyph'].text_font_size = '12pt'

    # Sankey with Drug names ===========================================================================================
    sankey = draw_sankey_plot(df=df_alt_lvl_min, col_id=col_sid,
                              cols=["Project_TCGA_More", "Sen_Level_Simple", "Res_Level_Simple"],
                              source_name="Tumor Type", target_name="Best sensitivity level",
                              other_threshold=0, edge_threshold=1,
                              edge_color="Target", label="",
                              cols_colors=[PTM2colors, Sen2colors, Res2colors],
                              width=args.width, height=args.height, label_text_font_size=9, data_aspect=3)
    sankey = sankey.opts(hooks=[hook])

    text_left = 'Tumor type'
    text_middle = 'Sensitvity levels'
    text_right = 'Resistance levels'
    sankey = sankey * hv.Text(0, 520, text=text_left).opts(text_font_size="14pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(500, 520, text=text_middle).opts(text_font_size="14pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(1000, 520, text=text_right).opts(text_font_size="14pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(500, 550, text="PRISM").opts(text_font_size="18pt", text_font_style="bold",
                                                              text_font=value("Helvetica"))

    sankey_state = hv.renderer('bokeh').get_plot(sankey).state
    sankey_state.output_backend = 'svg'
    export_svgs(sankey_state, filename=args.output)
    print("-file saved at %s" % args.output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cln_table', type=str, help='Table of clinical data.',
                        default="../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
    parser.add_argument('--alt_table', type=str, help='Table of aggregated alterations.',
                        default="../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv")
    parser.add_argument('--sam_table', type=str, help='Table of sample selection.',
                        default="../../../results/combined_alterations/selection/selection_samples_prism.tsv")
    parser.add_argument('--min_count', type=int, help='Min count per tumor type.', default=10)
    parser.add_argument('--width', type=int, help='Width of the output in pixels.', default=950)
    parser.add_argument('--height', type=int, help='Height of the output in pixels.', default=1250)
    parser.add_argument('--output', type=str, help='Path to output plot.',
                    default="../../../results/combined_alterations/other_plots/sankey_levels_prism.svg")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
