# -*- coding: utf-8 -*-
"""
@created: May 11 2022
@modified: May 11 2022
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Combine tables from different modalities (mutations, CNVs, fusions) into a single table.
"""

import argparse
import numpy as np
import os
import pandas as pd

# functions ============================================================================================================

def convert_num_to_str(x):
    try:
        y = "%d" % int(x)
    except:
        try:
            y = "%f" % float(x)
            if y=="nan":
                y = x
        except:
            y = x

    return y


def main(args):
    df_bio = pd.read_table(args.input_bio)
    df_cln = pd.read_table(args.input_cln)
    df_cnv = pd.read_table(args.input_cnv)
    df_fus = pd.read_table(args.input_fus)
    df_mut = pd.read_table(args.input_mut)

    # preprocess cnvs
    df_cnv["Alteration_Category"] = "c"

    # preprocess fusions
    df_fus = df_fus.rename(columns={"Fusion": "Hugo_Symbol"})
    df_fus["Alteration_Category"] = "f"

    # preprocess mutations
    df_mut["Alteration_Category"] = "m"

    # select columns required for each modality before concatenating
    col_tsb = "Tumor_Sample_Barcode"
    col_nsb = "Matched_Norm_Sample_Barcode"
    cols_req_cnv = [col_tsb, col_nsb, "Hugo_Symbol", "Alteration", "Alteration_Category"]
    cols_req_fus = [col_tsb, "Hugo_Symbol", "Gene_1", "Gene_2", "Alteration_Category"]
    cols_req_mut = [col_tsb, col_nsb, "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
                    "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "Hugo_Symbol", "Variant_Classification", "HGVSp", "HGVSc",
                    "all_effects", "Exon_Number", "Alteration_Category"]

    df_cnv = df_cnv[cols_req_cnv].copy()
    df_fus = df_fus[cols_req_fus].copy()
    df_mut = df_mut[cols_req_mut].copy()

    # add DNA_P
    df_cln["Sample_Id_DNA_P"] = df_cln[["Sample_Id_DNA_T", "Sample_Id_DNA_N"]].fillna("NA").apply("_vs_".join, axis=1)
    df_cnv["Sample_Id_DNA_P"] = df_cnv[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)
    df_mut["Sample_Id_DNA_P"] = df_mut[[col_tsb, col_nsb]].fillna("NA").apply("_vs_".join, axis=1)

    # it is important to apply here the selection of samples from cln_in_design in order to have at most one
    # pair (DNA, RNA) for each subject. alterations from the samples attached to the same subject may then be combined
    # for annotation by civic
    samples_rna = df_cln["Sample_Id_RNA_T"].drop_duplicates().unique().tolist()
    samples_dna = df_cln["Sample_Id_DNA_P"].drop_duplicates().unique().tolist()
    df_fus = df_fus.loc[df_fus["Tumor_Sample_Barcode"].isin(samples_rna)]
    df_cnv = df_cnv.loc[df_cnv["Sample_Id_DNA_P"].isin(samples_dna)]
    df_mut = df_mut.loc[df_mut["Sample_Id_DNA_P"].isin(samples_dna)]

    # concatenate
    df_alt = pd.concat((df_cnv[cols_req_cnv], df_fus[cols_req_fus], df_mut[cols_req_mut]), axis=0)

    # add Subject_Id, Civic_Disease and Sample_Type
    cols_cln = ["Subject_Id", "Tumor_Sample_Barcode", "Sample_Type", "Civic_Disease"]
    df_cln_rna = df_cln.loc[~df_cln["Sample_Id_RNA_T"].isnull()]
    df_cln_rna = df_cln_rna.rename(columns={"Sample_Id_RNA_T": "Tumor_Sample_Barcode"})
    df_cln_dna = df_cln.loc[~df_cln["Sample_Id_DNA_T"].isnull()]
    df_cln_dna = df_cln_dna.rename(columns={"Sample_Id_DNA_T": "Tumor_Sample_Barcode"})
    df_cln_rna_dna = pd.concat((df_cln_dna[cols_cln], df_cln_rna[cols_cln])).drop_duplicates()

    df_alt = df_alt.merge(df_cln_rna_dna, how="left", on="Tumor_Sample_Barcode")

    # save
    df_alt.to_csv(args.output, sep="\t", index=False)
    print("-file saved at %s" % args.output)


# run ==================================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine tables of CNVs, fusions and mutations.")
    parser.add_argument('--input_bio', type=str, help='Path to cnv table.',
                        default="../MetaPRISM/data/prism/clinical/curated/bio_prism_in_design_curated.tsv")
    parser.add_argument('--input_cln', type=str, help='Path to cnv table.',
                        default="../MetaPRISM/data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
    parser.add_argument('--input_cnv', type=str, help='Path to cnv table.',
                        default="./somatic_cna_civic_preprocess/all_samples.tsv")
    parser.add_argument('--input_fus', type=str, help='Path to cnv table.',
                        default="./somatic_fus_civic_preprocess/all_samples.tsv")
    parser.add_argument('--input_mut', type=str, help='Path to cnv table.',
                        default="./somatic_maf_civic_preprocess/all_samples.maf")
    parser.add_argument('--output', type=str, help='Path to output table.')
    args = parser.parse_args()

    for arg in vars(args):
        print("%s: %s" % (arg, getattr(args, arg)))
    print("\n")

    main(args)
