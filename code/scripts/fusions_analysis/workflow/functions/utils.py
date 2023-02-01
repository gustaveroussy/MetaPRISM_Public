import numpy as np
import pandas as pd

def load_table_samples(algo, cohort):
    filepath = "../../../data/%s/rna/%s/sample_list.tsv" % (cohort, algo)
    df_sam = pd.read_table(filepath, sep="\t")
    df_sam["Algo"] = algo
    return df_sam.loc[:,["Sample_Id", "Algo"]]


def extract_commonly_analyzed_samples(df_sam_a, df_sam_b):
    algos_a = sorted(df_sam_a["Algo"].unique())
    algos_b = sorted(df_sam_b["Algo"].unique())

    df_sam_a_agg = df_sam_a.sort_values(by="Algo").groupby("Sample_Id").agg({"Algo": "|".join}).reset_index()
    df_sam_b_agg = df_sam_b.sort_values(by="Algo").groupby("Sample_Id").agg({"Algo": "|".join}).reset_index()

    sam_a_all = df_sam_a_agg.loc[df_sam_a_agg["Algo"]=="|".join(algos_a), "Sample_Id"].values
    sam_b_all = df_sam_b_agg.loc[df_sam_b_agg["Algo"]=="|".join(algos_b), "Sample_Id"].values

    return list(set(sam_a_all).intersection(set(sam_b_all)))


def add_aggregated_callers(df, col_id, cols_call, algos_to_prefix):
    dfs_algo = []
    for algo, prefix in algos_to_prefix.items():
        col_call_algo = prefix + "_" + "Call"
        if col_call_algo in cols_call:
            df_algo = df.loc[df[col_call_algo]==1][[col_id]]
            df_algo["Algo"] = algo
            dfs_algo.append(df_algo)
    df_algo = pd.concat(dfs_algo,axis=0)
    df_algo = df_algo.sort_values(by="Algo").groupby(col_id).agg({"Algo":"|".join}).reset_index()
    return df.merge(df_algo, how="left", on=col_id)


