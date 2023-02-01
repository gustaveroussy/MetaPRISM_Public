import os
import numpy as np
import pandas as pd

def load_met_cov(dir_data, dir_models, names_features, names_models, names_selections, pool_level="imp"):
    pool_levels_ok = ["imp", "rep"]
    if not pool_level in pool_levels_ok:
        raise ValueError("-the value of pool_level should be in %s" % pool_levels_ok)

    dfs_met = []
    dfs_cov = []

    for name_features in names_features:
        filename_meta = "%s/processed/covs.final.tsv.gz" % name_features
        filepath_meta = os.path.join(dir_data, filename_meta)

        if os.path.exists(filepath_meta):
            df_meta = pd.read_table(filepath_meta)
            cols_rm = ["Min", "Max", "Redundancy"]
            df_meta = df_meta[[x for x in df_meta if x not in cols_rm]]

            for name_model in names_models:
                for name_selection in names_selections:
                    filename_met = "%s/%s_%s/mets.pooled_ax_%s.tsv.gz" % \
                            (name_features, name_selection, name_model, pool_level)
                    filename_cov = "%s/%s_%s/covs.pooled_ax_%s.tsv.gz" % \
                            (name_features, name_selection, name_model, pool_level)

                    filepath_met = os.path.join(dir_models, filename_met)
                    filepath_cov = os.path.join(dir_models, filename_cov)

                    if os.path.exists(filepath_met):
                        df_met = pd.read_table(filepath_met)
                        df_cov = pd.read_table(filepath_cov)

                        # replace covariate encoding by real names
                        df_cov = df_cov.rename(columns={"Covariate": "Code"})
                        df_cov = df_cov.merge(df_meta, how="left", on="Code")

                        dfs_met.append(df_met)
                        dfs_cov.append(df_cov)

    return pd.concat(dfs_met), pd.concat(dfs_cov)



def load_cov_ori(dir_data, names_features):
    dfs_cov_ori = []

    for name_features in names_features:
        filename_cov_ori = "%s/original/covs.original.tsv.gz" % name_features
        filepath_cov_ori = os.path.join(dir_data, filename_cov_ori)
        if os.path.exists(filepath_cov_ori):
            df_cov_ori = pd.read_table(filepath_cov_ori)
            df_cov_ori["Features"] = name_features
            dfs_cov_ori.append(df_cov_ori)

    return pd.concat(dfs_cov_ori)
