# -*- coding: utf-8 -*-
"""
@created: 13/09/21
@modified: 13/09/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Add annotations from FusionAnnotator to full fusions table.
"""

import ast
import argparse
import numpy as np
import pandas as pd
import re

# functions ============================================================================================================

def extract_fusion_type(x):
    if len(x)==0:
        return np.nan
    else:
        fusion_type = [e for e in x if e.startswith("INTRA") or e.startswith("INTER")]
        if len(fusion_type)==0:
            return np.nan
        else:
            fusion_type = fusion_type[0]
            regex = r"^[A-Za-z0-9]+"
            match = re.search(regex, fusion_type)
            if match is not None:
                return match.group(0)
            else:
                return np.nan

def simplify_annots(x):
    annots = [e for e in x if not e.startswith("INTRA") and not e.startswith("INTER")]
    if len(annots) == 0:
        return np.nan
    else:
        regex = r"^[A-Za-z0-9\_\-]+"
        annots_simple = []
        for annot in annots:
            if annot.startswith('["') and annot.endswith('"]'):
                annot = annot[2:-2]
            match = re.search(regex, annot)
            annot_simple = match.group(0)
            annots_simple.append(annot_simple)
        return "|".join(annots_simple)

def main(args):
    df_ann = pd.read_csv(args.input_annots, sep="\t", header=0)
    df_ann = df_ann.reset_index()
    df_ann.columns = ["Fusion_Id", "annots"]
    df_ann["annots"] = df_ann["annots"].apply(ast.literal_eval)
    df_ann["Fusion_Type"] = df_ann["annots"].apply(extract_fusion_type)
    df_ann["Annotations"] = df_ann["annots"].apply(simplify_annots)

    # add annotations
    df_fus = pd.read_csv(args.input_fusions, sep="\t", low_memory=False)
    cols_ann = ["Fusion_Id", "Fusion_Type", "Annotations"]
    df_fus = df_fus.merge(df_ann[cols_ann], how="left", on="Fusion_Id")

    # save
    df_fus.to_csv(args.output, sep="\t", index=False)

# parameters ===========================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add annotations from FusionAnnotator to fusion table.')
    parser.add_argument('--input_fusions', type=str, help="Input fusions file",
                       default="../../../data/tcga_6_samples/rna/fusions/tcga_6_samples_aggregated_callers.tsv")
    parser.add_argument('--input_annots', type=str, help="Input annotations file",
                       default="../../../data/tcga_6_samples/rna/fusions/tcga_6_samples_aggregated_FusionAnnotator_2.tsv")
    parser.add_argument('--output', type=str, help="Output file")
    args = parser.parse_args()

    print(args)
    main(args)
