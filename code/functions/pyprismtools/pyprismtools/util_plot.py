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

Util functions for defining plot parameters.
"""

import matplotlib
from matplotlib import cm
import re

def get_dict_labels_to_colors(labels_unique, sort=False, cmap_name=None, palette=None, alpha=None, colors_avoid=[]):
    labels2colors = dict()

    if palette is None and cmap_name is None:
        raise ValueError("Specify either palette or cmap_name")
    elif palette is None:
        cmap = cm.get_cmap(cmap_name)
        palette = [cmap(i) for i in range(cmap.N)]

    if sort:
        labels_unique = sorted(labels_unique)

    i = 0
    for label in labels_unique:
        color = palette[i % len(palette)]
        while color in colors_avoid:
            i = i+1
            color = palette[i % len(palette)]
        color = matplotlib.colors.to_hex(color)
        if alpha is not None:
            op = re.sub("0x", "", hex(round(alpha*255)))
            if len(op)==1:
                op = "0" + op
            color = color + op
        labels2colors[label] = color
        i = i+1

    return labels2colors
