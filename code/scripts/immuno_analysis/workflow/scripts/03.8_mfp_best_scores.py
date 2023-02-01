# -*- coding: utf-8 -*-
"""
@created: 24/06/22
@modified: 24/06/22
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Aggregated model quality scores from training on bagaev_2021 and tcga data.
"""

import argparse
import os
import numpy as np
import pandas as pd

# functions ============================================================================================================

def load_scores_model(folder, model):
    filepath = os.path.join(folder, model, "best/best_cv_scores.tsv")
    df_scores = pd.read_table(filepath)
    cols_scores = [x for x in df_scores if x not in ["type", "split"]]
    df_scores = df_scores.groupby("type").agg({x: "mean" for x in cols_scores})
    df_scores = df_scores.reset_index()
    df_scores.insert(0, "model", model)
    return df_scores

def main(args):
    models = [x for x in os.listdir(args.results_folder) if os.path.isdir(os.path.join(args.results_folder, x))]

    dfs_scores = []
    for model in models:
        df_scores = load_scores_model(folder=args.results_folder, model=model)
        dfs_scores.append(df_scores)
    df_scores = pd.concat(dfs_scores)

    # save
    df_scores.to_csv(args.output, sep="\t", index=False, float_format="%.3f")
    print("-file saved at %s" % args.output)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform prediction of Bagaev subtypes using pretrained model.')
    parser.add_argument('--results_folder', type=str,
                        help='Folder with CV results for each model trained.',
                        default="../../../results/immuno_analysis/tcga/mfp_models")
    parser.add_argument('--output', type=str, help="Path to output table",
                        default="../../../results/immuno_analysis/tcga/mfp_models/summary_best_cv_scores.tsv")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
