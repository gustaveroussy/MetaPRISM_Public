# -*- coding: utf-8 -*-
"""
@modified: Jul 27 2021
@created: Jul 27 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Useful functions for fitting models.
"""

import numpy as np
import os
import pandas as pd
import yaml
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, VotingClassifier
from sklearn.svm import SVC

def get_estimator(model):
    if model.startswith("MetaA"):
        return globals()[model.split("MetaA")[1]](estimators=[])
    elif model.startswith("MetaB"):
        return globals()[model.split("MetaB")[1]](estimators=[])
    else:
        if model=="LogisticRegressionStandard" or model=="LogisticRegressionLasso":
            return LogisticRegression()
        else:
            return globals()[model]()


def _proc_param_value(param_value):
    if param_value["func"]=="np.logspace":
        return np.logspace(start=param_value["start"], stop=param_value["stop"], num=param_value["num"])
    elif param_value["func"]=="np.linspace":
        return np.linspace(start=param_value["start"], stop=param_value["stop"], num=param_value["num"])
    elif param_value["func"]=="np.arange":
        return np.arange(start=param_value["start"], stop=param_value["stop"], step=param_value["step"])
    else:
        raise ValueError("Unknown value %s for 'func'" % param_value["func"])


def _get_param_grid(param_grid):
    for param_key, param_value in param_grid.items():
        if type(param_value)==dict:
            param_grid[param_key] = _proc_param_value(param_value)

    return param_grid


def _get_estimator_meta(folder, model):
    params_yaml = os.path.join(folder, model, "best/best_params.yaml")
    params = yaml.load(open(params_yaml), Loader=yaml.FullLoader)
    estimator = get_estimator(model)
    estimator = estimator.set_params(**params)

    return estimator


def _get_weight_meta(folder, models, weight_name):
    weights = []
    for model in models:
        filepath = os.path.join(folder, model, "best/best_cv_scores.tsv")
        df_scores = pd.read_csv(filepath, sep="\t")
        weight = df_scores.loc[df_scores["type"]=="test", weight_name].mean()
        weights.append(weight)
    return [weights]


def _get_param_grid_meta(param_grid):
    folder = "results/classification"
    models = list(param_grid["estimators"].keys())

    for model in models:
        if param_grid["estimators"][model]=="best_params":
            param_grid["estimators"][model] = _get_estimator_meta(folder, model)
    param_grid["estimators"] = [[(model, param_grid["estimators"][model]) for model in models]]

    for param_key, param_value in param_grid.items():
        if param_key == "weights" and type(param_value)==str:
            param_grid[param_key] = _get_weight_meta(folder, models, param_value)

        elif type(param_value)==dict:
            param_grid[param_key] = _proc_param_value(param_value)

    return param_grid


def get_param_grid(model, grids_yaml):
    param_grids = yaml.load(open(grids_yaml), Loader=yaml.FullLoader)

    if model not in param_grids.keys():
        raise ValueError("The model %s is not in the yaml file %s." % (model, grids_yaml))
    else:
        param_grid = param_grids[model]

    if model.startswith("MetaA") or model.startswith("MetaB"):
        return _get_param_grid_meta(param_grid)
    else:
        return _get_param_grid(param_grid)
