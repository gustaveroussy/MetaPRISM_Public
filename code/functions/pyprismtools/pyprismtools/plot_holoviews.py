# -*- coding: utf-8 -*-
"""
@created: Jun 24 2021
@modified: Jun 24 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions for drawing plot using the holoviews (https://holoviews.org/) library.
"""

import holoviews as hv
import matplotlib
import pandas as pd
import re

def prepare_sankey(df, col_id, cols, source_name="Source", target_name="Target", other_threshold=1, edge_threshold=1):
    df_c = df.copy()

    for col in cols:
        df_col_count = df_c[col].value_counts().sort_values()
        df_col_count = df_col_count.loc[df_col_count < other_threshold]
        index2other = {i: "Other" for i in df_col_count.index}
        df_c[col] = df_c[col].replace(index2other)

    df_sankeys = []

    for source, target in zip(cols, cols[1:]):
        df_nodup = df_c[[source, target] + [col_id]].drop_duplicates()
        df_nodup = df_nodup.fillna("Unknown")
        df_sankey = df_nodup.groupby([source, target])[col_id].nunique()
        df_sankey = df_sankey.reset_index(drop=False)
        df_sankey.columns = [source_name, target_name, "Value"]
        df_sankeys.append(df_sankey)

    # aggregate
    df_sankey = pd.concat(df_sankeys,axis=0)

    # remove cycles
    df_sankey = df_sankey.loc[df_sankey[source_name]!=df_sankey[target_name]]

    # remove edges with too low values
    df_sankey = df_sankey.loc[(df_sankey["Value"] >= edge_threshold)]

    return df_sankey


def draw_sankey_plot(df, col_id, cols, source_name="Source", target_name="Target", other_threshold=5, edge_threshold=5,
                     cols_colors=None, cmap="Set1", node_color="Source", edge_color="Target", label="",
                     label_position="outer", width=800, height=400, data_aspect=None, label_text_font_size=15):
    # prepare dataframe
    df_sankey = prepare_sankey(df, col_id, cols, source_name, target_name, other_threshold, edge_threshold)

    # prepare colors
    if cols_colors is not None:
        cmap = {}
        for colors in cols_colors:
            cmap.update(colors)

    node_color = source_name if node_color=="Source" else (target_name if node_color=="Target" else node_color)
    edge_color = source_name if edge_color=="Source" else (target_name if edge_color=="Target" else edge_color)

    sankey = hv.Sankey(df_sankey, label=label)
    sankey.opts(label_position=label_position, node_color=node_color, edge_color=edge_color, cmap=cmap)
    sankey.opts(width=width, height=height, data_aspect=data_aspect, label_text_font_size="%dpt" % label_text_font_size)
    return sankey
