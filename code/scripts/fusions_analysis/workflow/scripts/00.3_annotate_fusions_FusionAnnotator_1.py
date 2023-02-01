# -*- coding: utf-8 -*-
"""
@created: 13/09/21
@modified: 13/09/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Extract Fusion_Id column.
"""

import argparse
import pandas as pd

def main(args):
    df = pd.read_csv(args.input, sep="\t")
    df[["Fusion_Id"]].drop_duplicates().to_csv(args.output, index=False, sep="\t")

# parameters ===========================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract Fusion_Id column.')
    parser.add_argument('--input', type=str, help="Input file")
    parser.add_argument('--output', type=str, help="Output file")
    args = parser.parse_args()

    print(args)
    main(args)
