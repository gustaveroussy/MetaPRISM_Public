# -*- coding: utf-8 -*-
"""
@created: Aug 19 2020
@modified: Jul 06 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Select subset of a PanCancer MC3 MAF files based on a provided list of identifiers.
"""

import os
import pandas as pd
from   typing import Callable

DataFrame = pd.core.frame.DataFrame

def subset_maf_identifier(filepath_maf: str, col_identifier: str, identifiers: list, mod_identifier: Callable=None) -> DataFrame:
    df_maf_identifier = pd.read_csv(
        filepath_or_buffer = filepath_maf,
        sep                = "\t",
        usecols            = [col_identifier],
    )

    if mod_identifier is not None:
        df_maf_identifier.insert(0, "identifier", df_maf_identifier[col_identifier].apply(mod_identifier))
    else:
        df_maf_identifier = df_maf_identifier.rename(columns = {col_identifier: "identifier"})

    df_identifiers = pd.DataFrame({"identifier_specified": identifiers})
    df_maf_identifier = df_maf_identifier.merge(
        right    = df_identifiers,
        how      = "left",
        left_on  = "identifier",
        right_on = "identifier_specified"
    )
    indices_skip = df_maf_identifier.loc[df_maf_identifier["identifier_specified"].isnull()].index + 1

    df_maf = pd.read_csv(
        filepath_or_buffer = filepath_maf,
        sep                = "\t",
        skiprows           = indices_skip
    )

    return df_maf


def subset_maf_value(filepath_maf: str, col_value: str, value: str) -> DataFrame:
    df_maf_identifier = pd.read_csv(
        filepath_or_buffer = filepath_maf,
        sep                = "\t",
        usecols            = [col_value],
    )

    indices_skip = df_maf_identifier.loc[df_maf_identifier[col_value].apply(lambda x: value not in x)].index + 1

    df_maf = pd.read_csv(
        filepath_or_buffer = filepath_maf,
        sep                = "\t",
        skiprows           = indices_skip
    )

    return df_maf


def main(args):
    df_identifiers = pd.read_table(arg.identifiers)
    identifiers = df_identifiers["Aliquot_Id"].values.tolist()

    #### MC3 MAF CONTROLLED
    df_maf = subset_maf_identifier(
        filepath_maf = args.mc3_con,
        col_identifier = "Tumor_Sample_Barcode",
        identifiers    = identifiers,
        mod_identifier = lambda x: x
    )

    df_maf.to_csv(
        path_or_buf = args.mc3_con_sub,
        sep         = "\t",
        index       = False
    )

    print("-file saved at %s" %  args.mc3_con_sub)

    ## MC3 MAF PUBLIC

    df_maf = subset_maf_identifier(
        filepath_maf = args.mc3_pub,
        col_identifier = "Tumor_Sample_Barcode",
        identifiers    = identifiers,
        mod_identifier = lambda x: x
    )

    df_maf.to_csv(
        path_or_buf = args.mc3_pub_sub,
        sep         = "\t",
        index       = False
    )

    print("-file saved at %s" %  args.mc3_pub_sub)


# run ==================================================================================================================

if __name__ == "__main__":
    default_folder = "/Volumes/Data_II/Documents/data/TCGA/DNA/GRCH37/MAF/PAN-CANCER/MC3/files"

    parser = argparse.ArgumentParser(description="Perform a detail analysis of KRAS mutations in PAAD of TCGA.")
    parser.add_argument('--identifiers', type=str, help='',
                        default="../../../results/somatic_mutations/kras_checks/sam_bio_paad.tsv")
    parser.add_argument('--mc3_con', type=str, help='Path to MC3 controlled table.',
                        default="%s/mc3.v0.2.8.CONTROLLED.maf.gz" % default_folder)
    parser.add_argument('--mc3_pub', type=str, help='Path to MC3 public table.',
                        default="%s/mc3.v0.2.8.PUBLIC.maf.gz" % default_folder)
    parser.add_argument('--mc3_con_sub', type=str, help='Path to subsetted MC3 controlled table.',
                        default="../../../results/somatic_mutations/kras_checks/mc3.v0.2.8.CONTROLLED_PAAD.maf.gz")
    parser.add_argument('--mc3_pub_sub', type=str, help='Path to subsetted MC3 public table.',
                        default="../../../results/somatic_mutations/kras_checks/mc3.v0.2.8.PUBLIC_PAAD.maf.gz")
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n", end="")

    main(args)
