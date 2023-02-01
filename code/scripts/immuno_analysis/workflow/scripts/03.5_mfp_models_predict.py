# -*- coding: utf-8 -*-
"""
@created: 28/07/21
@modified: 02/08/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Uses a pretrained model to predict immune subtypes on new data.
"""

import argparse
import joblib
import numpy as np
import pandas as pd
import sys

from pyprism.data import load_table

sys.path.append("workflow/functions")
from getter import get_estimator, get_param_grid
from gridsearcher import GridSearchCVExtended
from transformer import ScalerPandas, mean_abs_dev
from utils import add_level_index, select_samples

# functions ============================================================================================================

def load_X(input_signatures):
    df_X = load_table(input_signatures, index_col=0)
    X = df_X.T
    return X


def main(args):
    X = load_X(args.input_signatures)
    X = add_level_index(X, args.input_annotation)
    X = select_samples(X, args.input_samples)

    pipe = joblib.load(args.input_model)
    Y_label = pipe.predict(X)
    Y_proba = pipe.predict_proba(X)
    classes = pipe.steps[-1][1].grid.best_estimator_.classes_.tolist()

    df_Y_label = pd.DataFrame(Y_label, index=X.index.get_level_values(0), columns=["Label"])
    df_Y_proba = pd.DataFrame(Y_proba, index=X.index.get_level_values(0), columns=["Proba_%s" % c for c in classes])
    df_Y = pd.concat((df_Y_label, df_Y_proba), axis=1)

    # save
    df_Y.to_csv(args.output, index=True, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform prediction of Bagaev subtypes using pretrained model.')
    parser.add_argument('--input_samples', type=str, nargs="*",
                        default="workflow/results/inputs/prism/selection_samples.tsv",
                        help='Table of samples to be retained.')
    parser.add_argument('--input_signatures', type=str,
                        default="workflow/results/inputs/prism/signatures_bagaev_2021_mfp.tsv",
                        help='Table of signatures scores per sample.')
    parser.add_argument('--input_annotation', type=str, default="pyprism_tcga",
                        help="Table of samples annotations containing 'Project_TCGA_More' column for preprocessing.")
    parser.add_argument('--input_model', type=str,
                        default="workflow/results/mfp_models/tcga/LogisticRegression/pipeline_train.joblib",
                        help="Path to the trained pipeline stored in a .joblib file.")
    parser.add_argument('--output', type=str,
                        default="workflow/results/inputs/tcga/mfp_subtypes_predicted_LogisticRegression.tsv",
                        help='Path to where trained pipeline will be saved.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
