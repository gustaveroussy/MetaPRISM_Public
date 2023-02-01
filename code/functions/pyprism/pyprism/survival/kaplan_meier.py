# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 21:24 2019

@author: Yoann Pradat

The present file implements functions for visualizing Kaplan-Meier curves
with different settings using matplotlib library.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dataclasses import dataclass
from sksurv.nonparametric import kaplan_meier_estimator
from typing import Tuple, List, Union

from .tests_survival import logrank_test

# type aliases
MplAxes = matplotlib.axes.Axes

@dataclass
class PlotKMConfig:
    """
    Settings for the plot_km function.

    Attributes
    ----------
    popu_label: str
        Name of original population
    mask_label: str
        Name of subpopulation
    time_col: str
        Name of the time-to-event column. Default:'time'
    evnt_col: str
        Name of the event column. Default:'event'
    use_popu_label: boolean
        If True, uses popu_label to label the curve. Default:False
    use_mask_label: boolean
        If True, uses mask_label to label the curve. Default:False
    linewidth: integer
        Linewidth of the curve. Default:2
    linestyle: str
        Linestyle of the curve. Default:'-'
    color: str
        Color of the line. Default:'orange'
    scatter_size: integer
        Size of the scatter points for censored individuals. Default:300
    """
    popu_label: str
    mask_label: str
    time_col: str='os'
    evnt_col: str='os_status'
    use_popu_label: bool=False
    use_mask_label: bool=False
    linewidth: int = 2
    linestyle: str='-'
    color: str = 'orange'
    scatter_size: int = 300

def plot_km(dataframe: pd.core.frame.DataFrame, config: PlotKMConfig, ax: MplAxes) -> bool:
    """
    This function plots the Kaplan-Meier estimate of the survival function for the input dataframe
    in the matplotlib Axes provided. Relies on kaplan_meier_estimator functions from sksurv.nonparametric

    Parameters
    ----------
    dataframe: pd.core.frame.DataFrame
        Data in the form of a dataframe with config.time_col and config.evnt_col column names
    config: PlotKMConfig
        Settings for the configuration of the plot.
    ax: matplotlib.axes._subplots.AxesSubplot
        Matplotlib axes in which the figure is drawn.

    Returns
    -------
    boolean stating whether the curves has been drawn or not.
    """
    if config.time_col not in dataframe.columns:
        raise ValueError("dataframe does not contain %s" % config.time_col)
    if config.evnt_col not in dataframe.columns:
        raise ValueError("dataframe does not contain %s" % config.evnt_col)

    df_nmos = dataframe[[config.evnt_col, config.time_col]].dropna()
    is_drawn = False

    if df_nmos.shape[0] == 0:
        print("""In populaton %s with mask %s, all individuals have one of
              missing '%s' or missing %s""" % (config.popu_label,config.mask_label, config.evnt_col, config.time_col))
    elif (df_nmos[config.evnt_col].nunique() == 1 and df_nmos[config.evnt_col].unique()[0] == False):
        print("""In populaton %s with mask %s, all individuals are
              censored """ % (config.popu_label,config.mask_label))
    elif df_nmos.shape[0] < 5:
        print("""In populaton %s with mask %s, there are less than 5
              individuals""" % (config.popu_label,config.mask_label))
    else:
        data_y = np.array([(event, time) for event, time in zip(df_nmos[config.evnt_col],
                                                                df_nmos[config.time_col])],
                          dtype=[('event', '?'), ('time', float)])

        time, prob = kaplan_meier_estimator(data_y['event'], data_y['time'])
        time = np.insert(time, 0, 0)
        prob = np.insert(prob, 0, 1.0)

        time_censored, prob_censored = [], []
        for t, p in zip(time, prob):
            if False in data_y[data_y['time'] == t]['event']:
                time_censored.append(t), prob_censored.append(p)

        if config.use_popu_label and config.use_mask_label:
            label = '%s %s %d/%d' % (config.popu_label,config.mask_label,data_y['event'].sum(),data_y.shape[0])
        elif config.use_popu_label:
            label = '%s %d/%d' % (config.popu_label,data_y['event'].sum(),data_y.shape[0])
        elif config.use_mask_label:
            label = '%s %d/%d' % (config.mask_label,data_y['event'].sum(),data_y.shape[0])
        else:
            label = '%d/%d' % (data_y['event'].sum(), data_y.shape[0])

        ax.scatter(time_censored, prob_censored, color=config.color, s=config.scatter_size, marker='+')
        ax.step(time, prob, color=config.color, linewidth=config.linewidth, linestyle=config.linestyle, label=label,
               where='post')
        is_drawn = True
    return is_drawn

@dataclass
class PlotCompareKMConfig:
    """
    subconfigs: List[PlotKMConfig]
        List of settings for each of the Kaplan-Meier estimates
    time_bounds: Tuple[float, float]
        Bounds on the x-axis of the time-to-event variable
    test_statistic: str
        Statistical test to use for comparing curves. For the moment only 'CMH_logrank' is available
        for more than 2 curves and 'linear_rank_logrank' for 2 curves. Default:'CMH_logrank'
    xy_statistic: Tuple[float, float]
        Coordinates in the axes fraction units to print the p-value of the test. Default: (0.6, 0.8)
    fontsize_statistic: Union[int,str]
        Fontsize to print the p-value of the test. Default:'small'
    """
    subconfigs: List[PlotKMConfig]
    time_bounds: Tuple[float, float]=(0, np.inf)
    test_statistic: str='CMH_logrank'
    xy_statistic: Tuple[float, float]=(0.6, 0.8)
    fontsize_statistic: Union[int, str]=15

def plot_compare_km(dataframes: List[pd.core.frame.DataFrame], config: PlotCompareKMConfig, ax: MplAxes) -> None:
    """
    This function plots the Kaplan-Meier estimates of the survival functions for each of the input dataframe
    in the matplotlib Axes provided. Relies on kaplan_meier_estimator functions from sksurv.nonparametric

    Parameters
    ----------
    dataframes: List[pd.core.frame.DataFrame]
        Data in the form of a dataframe with config.time_col and config.evnt_col column names
    config: PlotCompareKMConfig
        Settings for the configuration of the plot.
    ax: matplotlib.axes._subplots.AxesSubplot
        Matplotlib axes in which the figure is drawn.

    Returns
    -------
    None
    """
    are_drawn = []

    for df, subconfig in zip(dataframes, config.subconfigs):
        # drop NAs
        df  = df.dropna(axis=0, how='any', subset=[subconfig.time_col, subconfig.evnt_col])

        if df.shape[0] == 0:
            print('Population %s with mask %s has 0 individual' % (subconfig.popu_label, subconfig.mask_label))
        else:
            # select time bounds
            bounds_mask = ((df[subconfig.time_col] >= config.time_bounds[0]) &
                           (df[subconfig.time_col] <= config.time_bounds[1]))
            df_slct = df[[subconfig.evnt_col, subconfig.time_col]].dropna().loc[bounds_mask]
            # use plot function
            is_drawn=plot_km(df_slct, subconfig, ax)
            are_drawn.append(is_drawn)

    dfs_statistic = []
    evnt_cols = []
    time_cols = []
    n_groups = 0
    for df, is_drawn, subconfig in zip(dataframes, are_drawn, config.subconfigs):
        if is_drawn and df.shape[0] > 5 :
            dfs_statistic.append(df)
            evnt_cols.append(subconfig.evnt_col)
            time_cols.append(subconfig.time_col)
            n_groups += 1

    if n_groups >= 2:
        if config.test_statistic == 'CMH_logrank':
            statistic, df, p_value = logrank_test(dfs_statistic, evnt_cols, time_cols, type='CMH_logrank')
        elif config.test_statistic == 'Linear_rank_logrank':
            if n_groups > 2:
                raise ValueError('For test_statistic equal to Linear_rank_logrank only 2 groups are permitted')
            else:
                statistic, df, p_value = logrank_test(dfs_statistic, evnt_cols, time_cols, type='Linear_rank_logrank')
        elif config.test_statistic != None:
            raise ValueError("""Unsupported value of config.test_statistic \n Must be one of: None, CMH_logrank,
                             Linear_rank_logrank.""")

        if config.test_statistic != None:
            if p_value < 0.001:
                text = r'$\chi^2 P<0.001$'
            else:
                text = r'$\chi^2 P=%.3f$' % p_value
            ax.annotate(text=text, xy=config.xy_statistic, xycoords='axes fraction', fontsize=config.fontsize_statistic)
