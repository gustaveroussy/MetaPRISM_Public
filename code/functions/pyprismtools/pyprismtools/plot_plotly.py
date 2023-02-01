# -*- coding: utf-8 -*-
"""
@created: Nov 8 2021
@modified: Jun 02 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Util functions for drawing plotly plots with predefined themes.
"""

import pandas as pd
import plotly.graph_objects as go


def draw_barplot_counts(dfs, x, y, names2colors, textposition="outside", textfont=dict(size=10), textangle=90,
                        title=None, yaxis_range=None, barmode="group", sort=False, return_traces=False, name_upper=True,
                        **kwargs):
    if y != "Count":
        hovertemplate="<br>".join([
                "%{x}",
                y + ": " + "%{y:.3f}",
                "Count: " + "%{customdata}",
            ])
    else:
        hovertemplate="<br>".join([
                "%{x}",
                y + ": " + "%{y:}",
        ])

    if sort:
        items =  sorted(dfs.items())[::-1]
    else:
        items = dfs.items()

    data = []

    for name, df in items:
        df = df.sort_values(by=[x])
        if "Text" not in df.columns:
            text = df["Count"]
        else:
            text = df["Text"]

        if name_upper:
            name_bar = name.upper()
        else:
            name_bar = name

        bar = go.Bar(
            name=name_bar,
            x=df[x],
            y=df[y],
            text=text,
            textangle=textangle,
            textposition=textposition,
            textfont=textfont,
            customdata=df["Count"],
            hovertemplate=hovertemplate,
            marker_color=names2colors[name],
            **kwargs
        )
        data.append(bar)

    if return_traces:
        return data
    else:
        fig = go.Figure(data=data)

        fig = go.Figure(data=data)
        if yaxis_range is not None:
            fig.update_layout(yaxis_range=yaxis_range)
        if title is not None:
            fig.update_layout(barmode=barmode, title=dict(text=title), uniformtext_minsize=8, uniformtext_mode='show')
        return fig


def draw_barplot_counts_group_stack(dfs, names_stacks2colors, textfont_0=dict(color="black", size=10),
                                    textfont_1=dict(color="white", size=10), textposition_0="inside",
                                    textposition_1="outside", textangle=90, title=None, yaxis_range=None,
                                    return_traces=False, one_legend_per_stack=True, name_upper=True, **kwargs):
    hovertemplate="<br>".join([
        "%{x}",
        "%{y}",
    ])

    data = []
    stacks_seen = []
    showlegend_in_kwargs = "showlegend" in kwargs

    for i, (name, dt) in enumerate(dfs.items()):
        bases = None
        for stack, df in sorted(dt.items()):
            if bases is None:
                bases = [pd.Series(0, index=df.index)]
            bases.append(bases[-1] + df["Count"])

        for j, (stack, df) in enumerate(sorted(dt.items())[::-1]):
            if "Text" not in df:
                df["Text"] = df["Count"].astype(str)
                df["Text"] = df["Text"].replace("0","")

            if j==0:
                textfont = textfont_0
                textposition = textposition_0
            else:
                textfont = textfont_1
                textposition = textposition_1

            if one_legend_per_stack and stack not in stacks_seen:
                trace_name = stack
                showlegend = True
                stacks_seen.append(stack)
            else:
                if name_upper:
                    name_bar = name.upper()
                else:
                    name_bar = name
                trace_name = "%s - %s" % (name_bar, stack)
                showlegend = False

            if not showlegend_in_kwargs:
                kwargs["showlegend"] = showlegend

            bar = go.Bar(
                name=trace_name,
                x=df["Label"],
                y=df["Count"],
                text=df["Text"],
                textfont=textfont,
                textangle=textangle,
                textposition=textposition,
                hovertemplate=hovertemplate,
                offsetgroup=i,
                base=bases[-2-j],
                marker_color=names_stacks2colors[(name, stack)],
                **kwargs
            )

            data.append(bar)

    if return_traces:
        return data
    else:
        fig = go.Figure(data=data)
        if yaxis_range is not None:
            fig.update_layout(yaxis_range=yaxis_range)
        if title is not None:
            fig.update_layout(title=dict(text=title))
        return fig

