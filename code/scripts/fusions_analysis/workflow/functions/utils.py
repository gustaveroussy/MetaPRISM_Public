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


def merge_on_rows(df_x, df_y, how_rows="inner", prefer_y_over_x=True, prefer_exceptions=[],
                  fillna_where_possible=False, filepath_discrepancies=None, name_x="X", name_y="Y",
                  remove_na_from_discrepancies=False, write_only_discrepant_cols=False):
    rows_x = set(df_x.index)
    rows_y = set(df_y.index)
    rows_common = rows_x.intersection(rows_y)
    rows_x_only = rows_x.difference(set(rows_common))
    rows_y_only = rows_y.difference(set(rows_common))
    df_xn = df_x.copy()
    df_yn = df_y.copy()
    n_common = len(rows_common)

    if list(df_xn.index.names) == [None] or list(df_yn.index.names) == [None]:
        df_xn.index.names = ["Internal_Index"]
        df_yn.index.names = ["Internal_Index"]

    cols_x = set(df_xn.columns)
    cols_y = set(df_yn.columns)
    cols_common = cols_x.intersection(cols_y)
    cols_x_only = cols_x.difference(set(cols_common))
    cols_y_only = cols_y.difference(set(cols_common))

    def _add_level_columns(df, name, value, pos=0):
        df_cols = df.columns.to_frame()
        df_cols.insert(pos, name, value)
        df.columns = pd.MultiIndex.from_frame(df_cols)
        return df

    cols_common_drop = []
    dfs_discrepant = {}
    for col_c in cols_common:
        df_x_c = df_xn.loc[rows_common, col_c].fillna("NaN")
        df_y_c = df_yn.loc[rows_common, col_c].fillna("NaN")
        mask_discrepant = df_x_c != df_y_c
        n_discrepant = sum(mask_discrepant)

        if n_discrepant != 0:
            if remove_na_from_discrepancies:
                mask_discrepant = mask_discrepant & ~((df_x_c=="NaN") | (df_y_c=="NaN"))

            if sum(mask_discrepant) > 0:
                print("-warning: %s is discrepant for %d/%d common individuals" % \
                      (col_c, sum(mask_discrepant), n_common), flush=True)
            cols_common_drop.append(col_c)

            df_discrepant_x = df_xn.loc[rows_common].loc[mask_discrepant]
            if write_only_discrepant_cols:
                df_discrepant_x = df_discrepant_x[[col_c]]
            else:
                df_discrepant_x = df_discrepant_x[[c for c in df_discrepant_x.columns if c!=col_c] + [col_c]]
            df_discrepant_x = _add_level_columns(df_discrepant_x, name_x, name_x)

            df_discrepant_y = df_yn.loc[rows_common].loc[mask_discrepant]
            if write_only_discrepant_cols:
                df_discrepant_y = df_discrepant_y[[col_c]]
            else:
                df_discrepant_y = df_discrepant_y[[col_c] + [c for c in df_discrepant_y.columns if c!=col_c]]
            df_discrepant_y = _add_level_columns(df_discrepant_y, name_y, name_y)
            df_discrepant = pd.concat((df_discrepant_x, df_discrepant_y), axis=1).sort_index()

            dfs_discrepant["%s" % str(col_c)] = df_discrepant

    if filepath_discrepancies is not None and any([df.shape[0] > 0 for df in dfs_discrepant.values()]):
        with pd.ExcelWriter(filepath_discrepancies) as writer:
            for name, df in sorted(dfs_discrepant.items()):
                if df.shape[0] > 0:
                    one_saved = True
                    save_df_to_data(df, writer, sheet_name=name[:31], index=True, header=True, verbose=False)
            print("-excel file with discrepancies saved at %s" % filepath_discrepancies, flush=True)

    for col_c in cols_common_drop:
        if  col_c in prefer_exceptions:
            choose_y = not prefer_y_over_x
        else:
            choose_y = prefer_y_over_x

        if choose_y:
            del df_xn[col_c]
        else:
            del df_yn[col_c]
        cols_common.remove(col_c)

    index_names = df_xn.index.names
    df_xn = df_xn.reset_index(drop=False)
    df_yn = df_yn.reset_index(drop=False)
    cols_common = list(cols_common) + list(index_names)
    df_xn.loc[:, cols_common] = df_xn[cols_common].fillna("NaN")
    df_yn.loc[:, cols_common] = df_yn[cols_common].fillna("NaN")
    df_c = df_xn.merge(df_yn, how=how_rows, on=cols_common)
    df_c = df_c.replace(to_replace="NaN", value=np.nan)
    df_c = df_c.set_index(index_names)

    if list(df_c.index.names) == ["Internal_Index"]:
        df_c.index.names = [None]

    # fill in the gaps where possible
    if fillna_where_possible:
        for col_c in cols_common_drop:
            rows_na_x = df_x.loc[df_x[col_c].isnull()].index
            rows_common_na_x = list(set(rows_common).intersection(set(rows_na_x)))
            df_c.loc[rows_common_na_x, col_c] = df_y.loc[rows_common_na_x, col_c]

            rows_na_y = df_y.loc[df_y[col_c].isnull()].index
            rows_common_na_y = list(set(rows_common).intersection(set(rows_na_y)))
            df_c.loc[rows_common_na_y, col_c] = df_x.loc[rows_common_na_y, col_c]

    if len(rows_x_only) > 0 and how_rows in ["outer", "left"]:
        for col_x in cols_x:
            df_c.loc[rows_x_only, col_x]  = df_x.loc[rows_x_only, col_x]

    if len(rows_y_only) > 0 and how_rows in ["outer", "right"]:
        for col_y in cols_y:
            df_c.loc[rows_y_only, col_y]  = df_y.loc[rows_y_only, col_y]

    return df_c

