# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 2019

@author: Yoann Pradat

One data class and one function to configure matplotlib axes.
"""

import matplotlib as mpl

from dataclasses import dataclass, field
from typing import Dict, Union

# type alias
MplAxes = mpl.axes.Axes

def default_field(obj):
    return field(default_factory=lambda: obj)

@dataclass
class AxConfig:
    xlabel: Dict[str, Union[int, float, str]] = default_field({
		'text'       : 'xlabel',
        'fontsize'   : 20,
        'fontweight' : 'bold'
    })
    ylabel: Dict[str, Union[int, float, str]] = default_field({
		'text'       : 'ylabel',
        'fontsize'   : 20,
        'fontweight' : 'bold'
    })
    title: Dict[str, Union[int, float, str]] = default_field({
		'text'       : 'title',
        'fontsize'   : 25,
        'fontweight' : 'bold'
    })
    legend: Dict[str, Union[int, float, str]] = default_field({
        'show'     : True,
        'loc'      : 'best',
        'fontsize' : 20,
        'frameon'  : True
    })
    spines: Dict[str, Union[int, float, str]] = default_field({
		'bottom' : True,
        'left'   : True,
        'top'    : False,
        'right'  : False
    })
    xticks_params: Dict[str, Union[int, float, str]] = default_field({
		'which'     : 'major',
        'labelsize' : 12,
        'length'    : 3
    })
    yticks_params: Dict[str, Union[int, float, str]] = default_field({
		'which'     : 'major',
        'labelsize' : 12,
        'length'    : 3
    })
    xticks_params_ignore: bool=False
    yticks_params_ignore: bool=False


def set_axconfig(config: AxConfig, ax: MplAxes) -> None:
    # spines
    for key in config.spines.keys():
        ax.spines[key].set_visible(config.spines[key])
    # tick parameters
    if not config.xticks_params_ignore:
        ax.tick_params(
            axis='x',
            which=config.xticks_params['which'],
            labelsize=config.xticks_params['labelsize'],
            length=config.xticks_params['length']
        )
    if not config.yticks_params_ignore:
        ax.tick_params(
            axis='y',
            which=config.yticks_params['which'],
            labelsize=config.yticks_params['labelsize'],
            length=config.yticks_params['length']
        )
    # labels
    ax.set_xlabel(config.xlabel['text'],fontweight=config.xlabel['fontweight'],fontsize=config.xlabel['fontsize'])
    ax.set_ylabel(config.ylabel['text'],fontweight=config.ylabel['fontweight'],fontsize=config.ylabel['fontsize'])
    # title
    ax.set_title(label=config.title['text'],fontweight=config.title['fontweight'],fontsize=config.title['fontsize'])
    # legend
    if config.legend['show']:
        ax.legend(loc=config.legend['loc'], fontsize=config.legend['fontsize'], frameon=config.legend['frameon'])
