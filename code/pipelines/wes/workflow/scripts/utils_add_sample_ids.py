# -*- coding: utf-8 -*-
"""
@created: Jan 05 2022
@modified: Oct 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Add tumor and normal samples ids as columns to input file.
"""

import argparse
import pandas as pd

# functions ============================================================================================================


def main(args):
    # read
    df = pd.read_table(args.input)

    # add fields
    if args.tlabel in df:
        print("-WARNING: %s already in the columns of %s. Values are not overwritten." % (args.tlabel, args.table))
    else:
        df[args.tlabel] = args.tsample

    if args.nlabel in df:
        print("-WARNING: %s already in the columns of %s. Values are not overwritten." % (args.nlabel, args.table))
    else:
        df[args.nlabel] = args.nsample

    # write
    df.to_csv(args.output, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add Tumor_Sample and Normal_Sample fields.")
    parser.add_argument('--input', type=str, help='Path to input table.')
    parser.add_argument('--tsample', type=str, help='Name of tumor sample.')
    parser.add_argument('--nsample', type=str, help='Name of normal sample.')
    parser.add_argument('--tlabel', type=str, help='Column name for tumor sample.',
                        default="Tumor_Sample_Barcode")
    parser.add_argument('--nlabel', type=str, help='Column name for normal sample.',
                        default="Matched_Norm_Sample_Barcode")
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
