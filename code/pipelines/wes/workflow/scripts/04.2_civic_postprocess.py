# -*- coding: utf-8 -*-
"""
@created: Jan 24 2022
@modified: Jan 24 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Clean civic annotations.
"""

import argparse
import numpy as np
import os
import pandas as pd

# functions ============================================================================================================

def main(args):
    # load table
    df_table = pd.read_table(args.input)

    # select
    if df_table.shape[0]>0 and "CIViC_Matching_Disease" in df_table:
        df_table = df_table.loc[~df_table["CIViC_Matching_Disease"].isnull()].copy()
    else:
        df_table = df_table.iloc[:0,:]

    # save
    df_table.to_csv(args.output, index=False, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Postprocess CIViC annotations.")
    parser.add_argument('--input', type=str, help='Path to input table.')
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
