# -*- coding: utf-8 -*-
"""
@created: Oct 15 2021
@modified: Dec 30 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Sankey plots showing the distribution of treatments (drug names and classes) across tumor types.
"""

import argparse
import pandas as pd
import sys

# local functions
sys.path.append("workflow/functions")
from load import load_cln_cohorts
from util import get_table_counts

# pyprism
from pyprism.data import load_colors, load_table, load_cln
from pyprism.util import explode_df

# holoviews with bokeh backend
import holoviews as hv
from bokeh.io import export_svgs
from bokeh.core.properties import value, field
hv.extension('bokeh')

# pyprismtools
from pyprismtools import get_dict_labels_to_colors, draw_sankey_plot
from pyprismtools import draw_barplot_counts

# plotly
import plotly.io as pio
pio.templates.default = "plotly_white"


def main(args):
    # load cln table
    df_cln = load_cln_cohorts([args.cohort])

    # counts per tumor type
    df_counts_tt = df_cln["Project_TCGA_More"].value_counts()

    # drop na in drugs before biopsy and remove MISC - Not_TCGA
    df_cln = df_cln.dropna(subset=["Drugs_Before_Biopsy"])
    # df_cln = df_cln.loc[~df_cln["Project_TCGA_More"].isin(["MISC - Not_TCGA"])].copy()

    # extra only relevant columns
    cols_drug = ["Subject_Id", "Project_TCGA_More", "Biopsy_Selected", "Histological_Type", "MSKCC_Oncotree",
                 "Drugs_Before_Biopsy"]
    df_drug = explode_df(df=df_cln[cols_drug].dropna(subset=["Drugs_Before_Biopsy"]),
                         cols=["Drugs_Before_Biopsy"], sep="|")
    df_drug = df_drug.rename(columns={"Drugs_Before_Biopsy": "Drug_Name"})
    df_drug = df_drug.loc[df_drug["Drug_Name"]!="None"]

    # load drug table
    df_drug_table = load_table(args.drug_table)
    dci_2_class = {r["DCI"]: r["Class_Lvl_1"].replace("_", " ") for _,r in df_drug_table.iterrows()}

    # map to class
    df_drug["Drug_Class"] = df_drug["Drug_Name"].map(dci_2_class)

    # Colors
    PTM2colors = load_colors(sheet="Project_TCGA_More")
    DRC2colors = load_colors(sheet="Class_Lvl_1")
    DRC2colors = {k.replace("_", " "): v for k,v in DRC2colors.items()}
    DRN_unique = set(df_drug["Drug_Name"])
    DRN2colors = get_dict_labels_to_colors(labels_unique=DRN_unique, sort=True, cmap_name="tab20")

    # Select tumor types with sufficiently many tumors
    df_drug_min = df_drug.loc[df_drug["Project_TCGA_More"].isin(df_counts_tt[df_counts_tt>=args.min_count].index)]
    df_drug_rare = df_drug.loc[df_drug["Project_TCGA_More"].apply(lambda x: "Not_TCGA" in x)].copy()

    # For drug_min, group Not_TCGA subtypes under "Rare subtypes - Not_TCGA"
    mask_rare = df_drug_min["Project_TCGA_More"].apply(lambda x: "Not_TCGA" in x)
    mask_unkn = df_drug_min["Project_TCGA_More"]=="Unknown_Primary"
    type_rare = "Rare subtypes"
    type_unkn = "Unknown primary"
    df_drug_min.loc[mask_rare, "Project_TCGA_More"] = type_rare
    df_drug_min.loc[mask_unkn, "Project_TCGA_More"] = type_unkn

    # For drug_rare, replace MISC by MSKCC Oncotree
    mask_misc = df_drug_rare["Project_TCGA_More"]=="MISC - Not_TCGA"
    df_drug_rare.loc[mask_misc, "Project_TCGA_More"] = df_drug_rare.loc[mask_misc, "MSKCC_Oncotree"]

    # Hook for customizing fonts
    def hook(plot, element):
        plot.handles['text_1_glyph'].text_font = value('Helvetica')
        plot.handles['text_1_glyph'].text_font_size = '6pt'
        plot.handles['text_2_glyph'].text_font = value('Helvetica')
        plot.handles['text_2_glyph'].text_font_size = '6pt'

    # Sankey with Drug names ===========================================================================================
    sankey = draw_sankey_plot(df=df_drug_min, col_id="Subject_Id",
                              cols=["Project_TCGA_More", "Drug_Name"], source_name="Tumor Type", target_name="Drug",
                              other_threshold=5, edge_threshold=5,
                              edge_color="Source", label="",
                              cols_colors=[PTM2colors, DRN2colors],
                              width=args.widths[0], height=args.heights[0], label_text_font_size=6, data_aspect=10)
    sankey = sankey.opts(hooks=[hook])

    text_left = 'Number of sample-treatment\n relationships per tumor type'
    text_right = 'Number of sample-treatment\n relationships per drug'
    sankey = sankey * hv.Text(-100, 510, text=text_left).opts(text_font_size="9pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(1000, 510, text=text_right).opts(text_font_size="9pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))

    sankey_state = hv.renderer('bokeh').get_plot(sankey).state
    sankey_state.output_backend = 'svg'
    export_svgs(sankey_state, filename=args.outputs[0])
    print("-file saved at %s" % args.outputs[0])

    # Sankey with Drug class ===========================================================================================
    # Hook for customizing fonts
    def hook(plot, element):
        plot.handles['text_1_glyph'].text_font = value('Helvetica')
        plot.handles['text_1_glyph'].text_font_size = '8pt'
        plot.handles['text_2_glyph'].text_font = value('Helvetica')
        plot.handles['text_2_glyph'].text_font_size = '8pt'


    sankey = draw_sankey_plot(df=df_drug_min, col_id="Subject_Id",
                              cols=["Project_TCGA_More", "Drug_Class"], source_name="Tumor Type", target_name="Drug Class",
                              other_threshold=0, edge_threshold=10,
                              edge_color="Source", label="",
                              cols_colors=[PTM2colors, DRC2colors],
                              width=args.widths[1], height=args.heights[1], label_text_font_size=8, data_aspect=6)
    sankey.opts(hooks=[hook])

    text_left = 'Number of sample-treatment\n relationships per tumor type'
    text_right = 'Number of sample-treatment\n relationships per drug family'
    sankey = sankey * hv.Text(-100, 520, text=text_left).opts(text_font_size="10pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(1100, 520, text=text_right).opts(text_font_size="10pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))

    sankey_state = hv.renderer('bokeh').get_plot(sankey).state
    sankey_state.output_backend = 'svg'

    export_svgs(sankey_state, filename=args.outputs[1])
    print("-file saved at %s" % args.outputs[1])

    # Sankey with Drug class Rare tumor types ==========================================================================
    # Hook for customizing fonts
    def hook(plot, element):
        plot.handles['text_1_glyph'].text_font = value('Helvetica')
        plot.handles['text_1_glyph'].text_font_size = '6pt'
        plot.handles['text_2_glyph'].text_font = value('Helvetica')
        plot.handles['text_2_glyph'].text_font_size = '6pt'


    MSKCC2colors = get_dict_labels_to_colors(labels_unique=df_drug_rare.loc[mask_misc, "MSKCC_Oncotree"].unique(),
                                             cmap_name="tab20")
    TT2colors = {**PTM2colors, **MSKCC2colors}

    sankey = draw_sankey_plot(df=df_drug_rare, col_id="Subject_Id",
                              cols=["Project_TCGA_More", "Drug_Class"], source_name="Tumor Type", target_name="Drug Class",
                              other_threshold=0, edge_threshold=1,
                              edge_color="Source", label="",
                              cols_colors=[TT2colors, DRC2colors],
                              width=args.widths[2], height=args.heights[2], label_text_font_size=6, data_aspect=7)
    sankey.opts(hooks=[hook])

    text_left = 'Number of sample-treatment\n relationships per tumor type'
    text_right = 'Number of sample-treatment\n relationships per drug family'
    sankey = sankey * hv.Text(-100, 520, text=text_left).opts(text_font_size="8pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))
    sankey = sankey * hv.Text(1100, 520, text=text_right).opts(text_font_size="8pt", text_font_style="normal",
                                                              text_font=value("Helvetica"))

    sankey_state = hv.renderer('bokeh').get_plot(sankey).state
    sankey_state.output_backend = 'svg'
    export_svgs(sankey_state, filename=args.outputs[2])
    print("-file saved at %s" % args.outputs[2])

    # Side barplots for platinum salts ==================================================================================
    df_plat = df_drug_min.loc[df_drug_min["Drug_Class"]=="Platinum salts"]
    plat2colors = {"CISPLATINE": "#FFE6D5", "OXALIPLATINE":"#FFB380" , "CARBOPLATINE": "#FF6500"}
    x_stack_order = ["CARBOPLATINE", "OXALIPLATINE", "CISPLATINE"]

    # parameters
    x_group = "Drug_Class"
    x_stack = "Drug_Name"

    # counts
    stacks = df_plat[x_stack].dropna().unique()
    dfs = {stack: df_plat.loc[df_plat[x_stack]==stack] for stack in stacks}
    dfs = {k: get_table_counts(df, x_group, add_pct=True) for k, df in dfs.items()}
    dfs = {x: dfs[x] for x in x_stack_order}

    fig = draw_barplot_counts(dfs, x=x_group, y="Count", names2colors=plat2colors, sort=False, textposition="outside",
                              textangle=0, textfont=dict(size=24, family="Helvetica", color="black"),
                              yaxis_range=None, title="", barmode="stack")
    fig.update_xaxes(title_text="", showticklabels=False)
    fig.update_yaxes(title_font=dict(size=24, family="Helvetica", color="black"),
                     tickfont=dict(size=24, family="Helvetica", color="black"), dtick=250,
                     title_text="Number of\npatients", range=[0, 900])
    fig.update_layout(font_family="Helvetica", font_color="black", legend=dict(font=dict(size=20)), showlegend=False)

    # compute width and height
    width = 230
    height = 500

    fig.write_image(args.outputs[3], width=width, height=height, engine="kaleido")
    print("-file saved at %s" % args.outputs[3])

    # parameters
    # select only tumors analyzed in mutational signature analysis and from top tumor types
    x_group = "Project_TCGA_More"
    x_stack = "Drug_Name"

    df_sig = pd.read_table("../../../results/mutational_signatures/counts_mutations/counts_mutations_sbs_96_min_mut_prism.tsv")
    samples_sig = df_sig.columns[1:].tolist()
    df_cln_dna = load_cln(study=args.cohort).dropna(subset=["Sample_Id_DNA_T"])
    df_cln_dna = df_cln_dna.loc[df_cln_dna[x_group]!="MISC - Not_TCGA"]
    df_cnt_dna = df_cln_dna[x_group].value_counts()
    samples_min = df_cln_dna.loc[df_cln_dna[x_group].isin(df_cnt_dna[df_cnt_dna >= 10].index)]["Sample_Id_DNA_T"]
    samples_keep = list(set(samples_sig).intersection(set(samples_min)))
    subjects_keep = df_cln_dna.loc[df_cln_dna["Sample_Id_DNA_T"].isin(samples_keep)]["Subject_Id"]
    df_cnt_dna_keep = df_cnt_dna.loc[df_cnt_dna[df_cnt_dna >= 10].index]

    df_plat = df_drug_min.loc[df_drug_min["Drug_Class"]=="Platinum salts"]
    df_plat = df_plat.sort_values(by=["Subject_Id", "Drug_Name"])
    df_plat = df_plat.groupby("Subject_Id").agg({"Drug_Name": " & ".join, "Drug_Class": "first", "Project_TCGA_More":
                                                 "first"}).reset_index()

    plat2colors = {"CARBOPLATINE": "#F05365", "OXALIPLATINE":"#2EC4B6", "CISPLATINE": "#E5B3FE",
                   "CARBOPLATINE & CISPLATINE": "#FABC2A", "CARBOPLATINE & OXALIPLATINE": "#CBF3F0",
                   "CISPLATINE & OXALIPLATINE": "#FFDAB9", "CARBOPLATINE & CISPLATINE & OXALIPLATINE": "#F2EDEB"}

    df_plat_keep = df_plat.loc[df_plat["Subject_Id"].isin(subjects_keep)]
    df_plat_keep = df_plat_keep.groupby("Subject_Id").agg({x_stack: " & ".join, x_group: "first"}).reset_index()

    # counts
    stacks = df_plat_keep[x_stack].dropna().unique()
    dfs = {stack: df_plat_keep.loc[df_plat_keep[x_stack]==stack] for stack in stacks}
    dfs = {k: get_table_counts(df, x_group, add_pct=True) for k, df in dfs.items()}

    # correct percentage to account for patients that have not received any platinum
    for k in dfs.keys():
        dfs[k]["Percentage"] = dfs[k]["Count"]/df_cnt_dna.loc[dfs[k][x_group]].values
        dfs[k]["Count"] = ""
    dfs = {x: dfs[x] for x in plat2colors.keys() if x in dfs}

    fig = draw_barplot_counts(dfs, x=x_group, y="Percentage", names2colors=plat2colors, sort=False, textposition="inside",
                              textangle=0, textfont=dict(size=1, family="Helvetica", color="black"),
                              yaxis_range=None, title="", barmode="stack")
    fig.update_xaxes(title_text="", showticklabels=False, categoryorder="array",
                     categoryarray=list(df_cnt_dna_keep.index), linecolor="black", linewidth=1, mirror=True)
    fig.update_yaxes(title_font=dict(size=28, family="Helvetica", color="black"),
                     tickfont=dict(size=18, family="Helvetica", color="black"),
                     title_text="Patients, %", linecolor="black", linewidth=1, mirror=True)
    fig.update_layout(font_family="Helvetica", font_color="black",
                      legend=dict(y=1, traceorder='normal', font=dict(size=9)))

    # compute width and height
    width = 650
    height = 300

    fig.write_image(args.outputs[4], width=width, height=height, engine="kaleido")
    print("-file saved at %s" % args.outputs[4])

    fig.write_image(args.outputs[5], width=width, height=height, engine="kaleido")
    print("-file saved at %s" % args.outputs[5])


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make donut plots showing tumor type distribution for each cohort')
    parser.add_argument('--cohort', type=str, help='Names of the cohort.', default="prism")
    parser.add_argument('--min_count', type=int, help='Min count per tumor type.', default=10)
    parser.add_argument('--drug_table', type=str, help='Path to the drug table.',
                        default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
    parser.add_argument('--widths', type=int, nargs="+", help='Widths of the outputs in pixels.',
                        default=[390, 500, 500])
    parser.add_argument('--heights', type=int, nargs="+", help='Heights of the outputs in pixels.',
                        default=[1275, 650, 925])
    parser.add_argument('--outputs', type=str, nargs="+", help='Paths to output plots.',
                    default=["../../../results/data_overview/treatments/sankeys/sankey_tumor_type_drug_name_prism.svg",
                             "../../../results/data_overview/treatments/sankeys/sankey_tumor_type_drug_class_prism.svg",
                             "../../../results/data_overview/treatments/sankeys/sankey_tumor_type_drug_class_prism_rare.svg",
                             "../../../results/data_overview/treatments/sidebarplot_platinum.pdf",
                             "../../../results/data_overview/treatments/stacked_barplot_platinum.pdf",
                             "../../../results/figures_paper/F2b_mid.svg"])
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
