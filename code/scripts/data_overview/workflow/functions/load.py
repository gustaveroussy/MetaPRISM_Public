from pyprism.data import load_cln, load_bio

def _shorten(x, size=33):
    if type(x)==float:
        return x
    else:
        if len(x) > size:
            return x[:(size-3)//2] + "..." + x[-(size-3)//2:]
        else:
            return x

def _simplify(x):
    if type(x)==float:
        return x
    else:
        return x.split(",")[0]


def load_cln_cohorts(cohorts, shorten_size=23):
    """
    Load clinical data (subjects data) for each of the input cohorts into a dict of dataframes.
    """
    dfs_cln = {cohort: load_cln(study=cohort, mode="in_design") for cohort in cohorts}

    # perform some cosmetics on certain values to avoid lengthy names in plot labels
    cols_shorten = {"Primary_Site"}
    cols_simplify = {"Histological_Type"}
    cols_fillna = {"Primary_Site", "Project_TCGA", "Project_TCGA_More", "Histological_Type", "Histological_Type_Simple"}

    for name, df in dfs_cln.items():
        cols_df = set(df.columns)
        for col in cols_fillna.intersection(cols_df):
            df[col] = df[col].fillna("N/A").astype(str)
        for col in cols_shorten.intersection(cols_df):
            df["%s_Short" % col] = df[col].apply(_shorten, size=shorten_size)
        for col in cols_simplify.intersection(cols_df):
            df["%s_Simple" % col] = df[col].apply(_simplify)

    if len(cohorts)==1:
        return dfs_cln[cohorts[0]]
    else:
        return dfs_cln


def load_bio_cohorts(cohorts, shorten_size=23):
    """
    Load biospecimen data (samples data) for each of the input cohorts into a dict of dataframes.
    """
    # load dataframes with samples data into a dict of dataframes
    dfs_bio = {cohort: load_bio(study=cohort, mode="in_design") for cohort in cohorts}

    # perform some cosmetics on certain values to avoid lengthy names in plot labels
    cols_shorten = {"Biopsy_Site"}
    cols_simplify = {"Histological_Type"}
    cols_fillna = {"Biopsy_Site"}

    for name, df in dfs_bio.items():
        cols_df = set(df.columns)
        for col in cols_fillna.intersection(cols_df):
            df[col] = df[col].fillna("N/A").astype(str)
        for col in cols_shorten.intersection(cols_df):
            df["%s_Short" % col] = df[col].apply(_shorten, size=shorten_size)
        for col in cols_simplify.intersection(cols_df):
            df["%s_Simple" % col] = df[col].apply(_simplify)

    if len(cohorts)==1:
        return dfs_bio[cohorts[0]]
    else:
        return dfs_bio
