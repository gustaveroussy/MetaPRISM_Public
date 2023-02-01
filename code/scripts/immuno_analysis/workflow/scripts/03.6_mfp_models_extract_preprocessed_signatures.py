# -*- coding: utf-8 -*-
"""
@created: 06/08/21
@modified: 06/08/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Extract from a pretrained pipeline the transform signature expression matrix that is fed to the classifier.
"""

import argparse
import joblib
import numpy as np
import pandas as pd
import sys
from sklearn.pipeline import Pipeline

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
    X_t = Pipeline(pipe.steps[:-1]).transform(X)

    # save
    X_t.T.to_csv(args.output, index=True, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform prediction of Bagaev subtypes using pretrained model.')
    parser.add_argument('--input_samples', type=str, nargs="*",
                        default="workflow/results/inputs/prism/selection_samples.tsv",
                        help='Table of samples to be retained.')
    parser.add_argument('--input_signatures', type=str,
                        default="workflow/results/inputs/prism/signatures_bagaev_2021_mfp.tsv",
                        help='Table of signatures scores per sample.')
    parser.add_argument('--input_annotation', type=str, default="pyprism_prism",
                        help="Table of samples annotations containing 'Project_TCGA_More' column for preprocessing.")
    parser.add_argument('--input_model', type=str,
                        default="workflow/results/mfp_models/tcga/LogisticRegression/pipeline_train.joblib",
                        help="Path to the trained pipeline stored in a .joblib file.")
    parser.add_argument('--output', type=str,
                        default="workflow/results/inputs/prism/signatures_mfp_model_preprocessed.tsv",
                        help='Path to where trained pipeline will be saved.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
