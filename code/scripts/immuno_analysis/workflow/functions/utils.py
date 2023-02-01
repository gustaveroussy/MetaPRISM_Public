from pyprism.data import load_table, load_cln

def add_level_index(X, input_annotation, col=["Project_TCGA_More", "TCGA_project"]):
    if input_annotation.startswith("pyprism_"):
        cohort = input_annotation.split("pyprism_")[1]
        df_Y = load_cln(study=cohort)
        df_Y = df_Y.dropna(subset=["Sample_Id_RNA_T"]).set_index("Sample_Id_RNA_T")
    else:
        df_Y = load_table(input_annotation, index_col=0)

    col_on = None
    for c in col:
        if c in df_Y:
            col_on = c
            break

    if col_on is None:
        raise ValueError("None of the col names specified was found in the table %s." % input_annotation)

    Z = df_Y[col_on]
    X = X.merge(Z, left_index=True, right_index=True, how="left")
    X = X.set_index(col_on, append=True)

    return X


def select_samples(X, input_samples):
    if input_samples is None or input_samples=="" or input_samples==[]:
        return X
    else:
        if type(input_samples)==list:
            S = load_table(input_samples[0])
        else:
            S = load_table(input_samples)
        S = S.loc[S["Use_mfp_model"]==1]
        return X.loc[S["Sample_Id"]]
