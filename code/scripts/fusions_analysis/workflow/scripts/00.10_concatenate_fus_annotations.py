# -*- coding: utf-8 -*-
"""
@created: 04 Mar 22
@modified: 04 Mar 22
@authors: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Aggregate mutations annotations from CIViC, OncoKB and CGI (to be done).
"""

import argparse
import gzip
import numpy as np
import pandas as pd
import re
import sys
sys.path.append("../pipeline_cln/workflow/functions")

from util_curate import merge_on_rows

# functions ============================================================================================================

def correct_dtype_from_num_to_str(x):
    try:
        y = "%d" % int(float(x))
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = np.nan
        except:
            y = x

    return y


def main(args):
    df_okb = pd.read_table(args.okb)
    df_civ = pd.read_table(args.civ)

    # reorder columns civ
    cols_ppd = [x for x in df_civ if re.search("^Prognostic:|^Predictive:|^Diagnostic:", x) is not None]
    cols_oth = [x for x in df_civ if x not in cols_ppd]
    df_civ = df_civ[cols_oth+cols_ppd]

    suffixes = ["Breakpoint_1", "Breakpoint_2", "Reads", "Chr_1", "Chr_2"]
    for df in [df_okb, df_civ]:
        cols = [x for x in df if any([x.endswith(s) for s in suffixes])]
        for col in cols:
            df[col] = df[col].apply(correct_dtype_from_num_to_str)

    # info messages
    okb_and_civ = len(set(df_okb["Row_Id"]).intersection(set(df_civ["Row_Id"])))
    okb_not_civ = len(set(df_okb["Row_Id"]).difference(set(df_civ["Row_Id"])))
    civ_not_okb = len(set(df_civ["Row_Id"]).difference(set(df_okb["Row_Id"])))
    okb_or_civ = len(set(df_civ["Row_Id"]).union(set(df_okb["Row_Id"])))
    print("-INFO: there are %d rows in OncoKB and CIViC" % okb_and_civ)
    print("-INFO: there are %d rows in OncoKB not in CIViC" % okb_not_civ)
    print("-INFO: there are %d rows in CIViC not in OncoKB" % civ_not_okb)
    print("-INFO: there are %d rows in OncoKB or CIViC" % okb_or_civ)

    df_ann = merge_on_rows(df_x=df_okb.set_index("Row_Id"),
                           df_y=df_civ.set_index("Row_Id"),
                           how_rows="outer",
                           prefer_y_over_x=True,
                           prefer_exceptions=[],
                           fillna_where_possible=True,
                           filepath_discrepancies=None,
                           name_x="OncoKB",
                           name_y="CIViC",
                           remove_na_from_discrepancies=True,
                           write_only_discrepant_cols=True)

    df_ann = df_ann.reset_index(drop=True)

    # check and save
    assert df_ann.shape[0] == okb_or_civ
    df_ann = df_ann.replace("nan", np.nan)
    df_ann.to_csv(args.output, index=False, sep="\t")

# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Aggregate mutations annotations.")
    parser.add_argument("--civ", type=str, help="Path to CIViC-annotated mutations table.",
                      default="../../../data/prism/rna/fusions/prism_annotated_filtered_civic.tsv.gz")
    parser.add_argument("--okb", type=str, help="Path to OncoKB-annotated mutations table.",
                      default="../../../data/prism/rna/fusions/prism_annotated_filtered_oncokb.tsv.gz")
    parser.add_argument('--output', type=str, help='Path to output table.',
                      default="../../../data/prism/rna/fusions/prism_annotated_filtered_union_ann.tsv.gz")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
