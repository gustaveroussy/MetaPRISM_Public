# -*- coding: utf-8 -*-
"""
@created: Oct 26 2022
@modified: Oct 26 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Remove empty fields from input table.
"""

import argparse
import pandas as pd

# functions ============================================================================================================


def main(args):
    # read
    df = pd.read_table(args.input)

    # remove empty fields
    empty_level = max(0, min(args.level, 1))
    all_fields = df.columns.tolist()
    empty_fields = df.columns[df.isnull().mean(axis=0)>=empty_level].tolist()
    keep_fields = [x for x in all_fields if x not in empty_fields]
    print("-INFO: identified %d/%d completely empty fields" % (len(empty_fields), len(all_fields)))
    if len(empty_fields)>0:
        print("\t" + "\n\t".join(empty_fields))
    df = df[keep_fields]

    # write
    df.to_csv(args.output, sep="\t", index=False)

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove empty fields.")
    parser.add_argument('--input', type=str, help='Path to input table.')
    parser.add_argument('--level', type=float, help='Level of emptiness above which fields are removed.', default=1)
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
