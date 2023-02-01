# -*- coding: utf-8 -*-
"""
@created: 27 Jul 2021
@modified: 18 Feb 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Learns a model for performing immune subtyping as presented in the paper PMID: 34019806.
"""

import argparse
import joblib
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import yaml

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, StratifiedKFold

from sklearn.metrics import accuracy_score, recall_score, precision_score, f1_score, roc_auc_score
from sklearn.metrics import make_scorer

from pyprism.data import load_table

sys.path.append("workflow/functions")
from getter import get_estimator, get_param_grid
from gridsearcher import GridSearchCVExtended
from transformer import SubpartsScalerPandas, mean_abs_dev
from utils import add_level_index

# functions ============================================================================================================

def load_XY(input_signatures, input_annotation):
    df_X = load_table(input_signatures, index_col=0)
    df_Y = load_table(input_annotation, index_col=0)
    X = df_X.T
    Y = df_Y.loc[X.index, "MFP"]

    return X,Y


def main(args):
    X, Y = load_XY(args.input_signatures, args.input_annotation)
    X = add_level_index(X, args.input_annotation, col=["Project_TCGA_More", "TCGA_project"])
    kfold_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=args.random_state)

    scoring = {"accuracy": make_scorer(accuracy_score),
               "recall": make_scorer(recall_score, average="macro", zero_division=0),
               "precision": make_scorer(precision_score, average="macro", zero_division=0),
               "f1_weighted": make_scorer(f1_score, average="weighted", zero_division=0),
               "f1_macro": make_scorer(f1_score, average="macro", zero_division=0)}

    pipe = Pipeline([
        ('scaler_median', SubpartsScalerPandas(center_func=np.median, scale_func=mean_abs_dev, clip=None)),
        ('gridsearch', GridSearchCVExtended(
            grid = GridSearchCV(
                estimator=get_estimator(args.model_name),
                param_grid=get_param_grid(args.model_name, args.grids_yaml, args.cohort),
                cv=kfold_cv,
                scoring=scoring,
                refit="f1_weighted",
                return_train_score=True,
                verbose=1,
                n_jobs=args.threads
            ),
            save_best_estimators_cv=True,
            save_best_estimator_train=False,
            save_best_cv_scores=True,
            save_best_params=True,
            save_folder=str(Path(args.output).parent.absolute()),
            verbose=True
        )),
    ])

    pipe = pipe.fit(X, Y)

    # save
    with open(args.output, "wb") as file:
        joblib.dump(pipe, file)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Learn a classification model using subtypes from Bagaev supp.')
    parser.add_argument('--input_signatures', type=str, default="resources/bagaev_2021/Pan_TCGA/signatures.tsv",
                        help='Table of signatures scores per sample.')
    parser.add_argument('--input_annotation', type=str, default="resources/bagaev_2021/Pan_TCGA/annotation.tsv",
                        help="""Table of samples annotations containing immune subtypes in 'MFP' column and
                        'Project_TCGA_More' column for preprocessing.""")
    parser.add_argument('--cohort', type=str, default="bagaev_2021",
                        help="""Name of the cohort. Useful for find best parameters from previous models when running
                        a meta model.""")
    parser.add_argument('--model_name', type=str, default="LogisticRegression", help='Class of model to be used.')
    parser.add_argument('--grids_yaml', type=str, default="config/mfp_models_grids.yaml",
                        help='Path to yaml file defining the param grids.')
    parser.add_argument('--random_state', type=int, default=123, help='Seed random generator.')
    parser.add_argument('--threads', type=int, default=4, help='Number of jobs available.')
    parser.add_argument('--output', type=str, help='Path to where trained pipeline will be saved.',
                default="../../../results/immuno_analysis/bagaev_2021/mfp_models/LogisticRegression/pipeline_train.joblib")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
