# -*- coding: utf-8 -*-
"""
@created: Nov 8 2021
@modified: Nov 8 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Functions util throughout the data analysis.
"""

import re
import pandas as pd

from pyprism.data import load_bio, load_cln

def preprocess_wes_mut(df_mut, cohort, cols_bio=[], cols_cln=[], select_pairs=False, selection_mut="all", verbose=True):
    """Add bio and/or cln attributes to the mutations table and perform a selection of one pair tumor/normal per subject
    if requested.

    Parameters
    ----------
    df_mut: dataframe,
        Table of mutations.
    cohort: str,
        Name of the cohort
    cols_bio: list, optional
        List of bio attributes.
    cols_cln: list, optional
        List of cln attributes.
    select_pairs: bool, optional
        If set to True, 1 pair tumor/normal is selected for each subject.
    selection_mut: str, optional
        Here are the available modes
            - 'all': no selection of variants
            - 'non_synonymous': select variants for which the value of 'Variant_Classification' is one of
                * Frame_Shift_Del
                * Frame_Shift_Ins,
                * Splice_Site
                * Translation_Start_Site
                * Nonsense_Mutation
                * Nonstop_Mutation
                * In_Frame_Del
                * In_Frame_Ins
                * Missense_Mutation
                * Start_Codon_Del
                * Start_Codon_SNP
                * Stop_Codon_Del
                * Stop_Codon_Ins
            - 'truncating': select variants for which the value of 'Variant_Classification' is one of
                * Frame_Shift_Del
                * Frame_Shift_Ins,
                * Splice_Site
                * Nonsense_Mutation
    verbose: bool, optional
        Should info messages be printed?

    Returns
    -------
    A dataframe of mutations with possibly additional attributes and possibly a reduced number of lines.
    """
    df_muc = df_mut.copy()

    # add Tumor_Sample_Id, Normal_Sample_Id
    if cohort=="tcga":
        df_muc["Tumor_Sample_Id"] = df_muc["Tumor_Sample_Barcode"].str[:20]
        df_muc["Normal_Sample_Id"] = df_muc["Matched_Norm_Sample_Barcode"].str[:20]
    else:
        df_muc["Tumor_Sample_Id"] = df_muc["Tumor_Sample_Barcode"]
        df_muc["Normal_Sample_Id"] = df_muc["Matched_Norm_Sample_Barcode"]

    # bio attributes
    cols_bio_min = ["Sample_Id", "Subject_Id"]
    cols_bio_use = list(set(cols_bio).union(set(cols_bio_min)))
    cols_bio_rmv = list(set(cols_bio_use).difference(set(cols_bio)))
    df_bio = load_bio(cohort)[cols_bio_use].drop_duplicates()
    df_muc = df_muc.merge(df_bio, how="left", left_on="Tumor_Sample_Id", right_on="Sample_Id")

    # cln attributes
    cols_cln_min = ["Subject_Id", "Sample_Id_DNA_T", "Sample_Id_DNA_N"]
    cols_cln_use = list(set(cols_cln).union(set(cols_cln_min)))
    cols_cln_rmv = list(set(cols_cln_use).difference(set(cols_cln)))
    df_cln = load_cln(cohort)[cols_cln_use].drop_duplicates()
    df_cln = df_cln.dropna(subset=cols_cln_min, how="any")
    df_muc = df_muc.merge(df_cln, how="left", on="Subject_Id")

    # tumor/normal pairs selection
    if select_pairs:
        col_tsid_m = "Tumor_Sample_Id"
        col_tsid_c = "Sample_Id_DNA_T"
        col_nsid_m = "Normal_Sample_Id"
        col_nsid_c = "Sample_Id_DNA_N"
        col_pid = "Tumor_Vs_Normal"

        df_muc[col_pid] = df_muc[[col_tsid_m, col_nsid_m]].apply("_vs_".join, axis=1)
        df_cln[col_pid] = df_cln[[col_tsid_c, col_nsid_c]].apply("_vs_".join, axis=1)
        df_muc = df_muc.loc[df_muc[col_pid].isin(df_cln[col_pid])]

        # info messages
        if verbose:
            print("-number of unique pairs: %d" % len(set(df_muc[col_pid])), flush=True)
            print("-number of unique subjects: %d" % len(set(df_muc["Subject_Id"])), flush=True)

        del df_muc[col_pid]

    # remove unwanted attributes
    cols_rmv = list((set(cols_bio_rmv).union(set(cols_cln_rmv))).difference(set(cols_bio).union(set(cols_cln))))
    cols_rmv = list(set(cols_rmv).union(set(["Tumor_Sample_Id", "Normal_Sample_Id"])))
    df_muc = df_muc[[x for x in df_muc if x not in cols_rmv]]

    # variants selection
    if selection_mut not in ["all", "annotated"]:
        if selection_mut=="non_synonymous":
            col_keep = "Variant_Classification"
            val_keep = ["Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site",
                        "Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                        "Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins","Translation_Start_Site"]
        elif selection_mut=="truncating":
            col_keep = "Variant_Classification"
            val_keep = ["Frame_Shift_Del","Frame_Shift_Ins","Splice_Site", "Nonsense_Mutation", "Splice_Site",
                        "Splice_Region", "Nonstop_Mutation", "Stop_Codon_Del", "Stop_Codon_Ins"]
        else:
            raise ValueError("-unsupported value %s for argument 'selection_mut'." % selection_mut)

        df_muc = df_muc.loc[df_muc[col_keep].isin(val_keep)]

        if verbose:
            print("-number of selected variants: %d" % df_muc.shape[0], flush=True)

    return df_muc



def split_targeted_therapy_targets(x, sep="|"):
    """
    Split multiple-targets therapies into multiple single-target therapies.

    Parameters
    ----------
    x: str,
        a string containing entries (e.g drug names) separated by sep
    sep: str, default="|"
        the separator

    Return
    ------
    str
    """
    classes = x.split(sep)
    prefix = "Targeted_Therapy - "
    classes_new = []
    classes_old = []
    for c in classes:
        if c.startswith(prefix):
            targets = c.split(prefix)[1].split("/")
            if len(targets) > 1:
                classes_new += ["%s%s" % (prefix, target) for target in targets]
                classes_old += [c]
    classes = list((set(classes).difference(set(classes_old))).union(set(classes_new)))
    return sep.join(classes)

