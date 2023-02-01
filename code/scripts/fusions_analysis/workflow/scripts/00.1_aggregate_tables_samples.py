# -*- coding: utf-8 -*-
"""
@created: 27/08/21
@modified: 22/12/21
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

Create aggregated tables of fusions and list of samples processed.
"""

import argparse
import os
import pandas as pd
from tqdm import tqdm

# function =============================================================================================================

def get_corrected_sample_id(sample_id):
    if sample_id=="M727-T1-ARN":
        return "M727RE-T1-ARN"
    elif sample_id=="MR012-T-ARN":
        return "MR012-T2-ARN"
    else:
        return sample_id


def get_sample_id_from_rna_fus(x, cohort=None):
    if cohort == "prism":
        return get_corrected_sample_id(x+"-ARN")
    elif cohort == "met500":
        return x
    elif cohort == "tcga_validation":
        return x
    else:
        raise ValueError("Unsupported value '%s' for 'cohort'. Choose 'met500', 'prism' or 'tcga_validation'." % cohort)


def main(args):
    folder = "%s/%s/rna/%s" % (args.data_folder, args.cohort, args.algo_folder)
    print("-processing algo %s" % args.algo_folder)

    if os.path.exists(args.output_list):
        df_sam = pd.read_table(args.output_list)
    else:
        subfolder_list = [x for x in os.listdir(folder) if os.path.isdir("%s/%s" % (folder, x))]
        df_sam = pd.DataFrame({"Sample_Id": subfolder_list})
        df_sam["Sample_Id"] = df_sam["Sample_Id"].apply(get_sample_id_from_rna_fus, cohort=args.cohort)
        df_sam.to_csv(args.output_list, sep="\t", index=False, header=False)
        df_sam["Isilon_Id"] = df_sam["Sample_Id"]

    isilon_list = df_sam["Isilon_Id"].tolist()
    sample_list = df_sam["Sample_Id"].tolist()

    df_list = []
    for isilon, sample in tqdm(zip(isilon_list, sample_list), total=len(isilon_list)):
        if os.path.exists("%s/%s" % (folder, isilon)):
            files = [x for x in os.listdir("%s/%s" % (folder, isilon)) if ".tsv" in x or ".txt" in x]
            files = [x for x in files if isilon in x]
            if len(files) > 0:
                file = files[0]
                df = pd.read_table("%s/%s/%s" % (folder, isilon, file), sep="\t", low_memory=False)
                df["Sample_Id"] = sample
                df_list.append(df)

    df_algo = pd.concat(df_list, axis=0)
    df_algo.to_csv(args.output_agg, index=False, sep="\t", compression="gzip")
    print("-saved aggregated table for algo %s" % args.algo_folder)


# parameters ===========================================================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Aggregate fusion calls across samples for a given caller.')
    parser.add_argument('--cohort', type=str, help="Cohort name", default="tcga_validation")
    parser.add_argument('--algo_folder', type=str, help="Name of folder where algo data is located.",
                        default="arriba")
    parser.add_argument('--data_folder', type=str, help="Path to data folder", default="../../../data")
    parser.add_argument('--output_list', type=str, help="Path to output table.",
                        default="../../../data/tcga_validation/rna/arriba/sample_list.tsv")
    parser.add_argument('--output_agg', type=str, help="Path to output table.",
                        default="../../../data/tcga_validation/rna/arriba/tcga_validation_arriba.tsv.gz")
    args = parser.parse_args()

    print(args)
    main(args)
