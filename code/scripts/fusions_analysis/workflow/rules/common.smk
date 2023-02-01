from itertools import product
import ntpath
import pandas as pd
from snakemake.utils import min_version
import yaml

min_version("5.4.0")
configfile: "config/config.yaml"

D_FOLDER = config["filepaths"]["base"]
F_FOLDER = "../../../results/figures_paper"
R_FOLDER = "../../../results/fusions_analysis"
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)

COHORTS = config["data"]["cohorts"]
COHORTS_NOT_TCGA = [x for x in COHORTS if x.upper()!="TCGA"]
COHORTS_NOT_VAL = [x for x in COHORTS if x.upper()!="TCGA_VALIDATION"]
ALGOS = list(set().union(*[config["data"]["aggregate"][cohort]["algos"] for cohort in COHORTS]))

##### Helper functions #####

def filter_combinator(combinator, comblist, white_list=True):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
        # Use frozenset instead of tuple
        # in order to accomodate
        # unpredictable wildcard order
            if white_list:
                if frozenset(wc_comb) in comblist:
                    yield wc_comb
            else:
                if frozenset(wc_comb) not in comblist:
                    yield wc_comb
    return filtered_combinator



def get_allowed_doublets_cohort_algo():
    allowed = []
    for cohort in COHORTS:
        for algo in config["data"]["aggregate"][cohort]["algos"]:
            allowed.append(frozenset({("cohort", cohort), ("algo", algo)}))
        return filter_combinator(product, allowed, white_list=True)


def get_allowed_algos(wildcards):
    algos = []
    return config["data"]["aggregate"][wildcards.cohort]["algos"]


def get_input_aggregate_tables_callers(wildcards):
    if wildcards.cohort in COHORTS_NOT_TCGA:
        return expand("%s/{cohort}/rna/{algo}/{cohort}_{algo}.tsv.gz" % D_FOLDER,
                      cohort=wildcards.cohort, algo=config["data"]["aggregate"][wildcards.cohort]["algos"])
    else:
        algos = config["data"]["aggregate"][wildcards.cohort]["algos"]
        files = [FILEPATHS[wildcards.cohort]["rna_fus"][algo] for algo in algos]
        return ["%s/%s" % (D_FOLDER, file) for file in files]


def get_table(cohort):
    filepath_cln = "%s/%s/clinical/curated/cln_%s_in_design_curated.tsv" % (D_FOLDER, cohort, cohort)
    filepath_bio = "%s/%s/clinical/curated/bio_%s_in_design_curated.tsv" % (D_FOLDER, cohort, cohort)
    df_cln = pd.read_table(filepath_cln)
    df_bio = pd.read_table(filepath_bio)
    cols_cln = ["Subject_Id", "Project_TCGA_More", "MSKCC_Oncotree", "Civic_Disease"]
    table = df_bio.merge(df_cln[cols_cln], how="left", on="Subject_Id")
    return table.set_index("Sample_Id")


def get_column_table_sample(wildcards, col):
    """Get the value of the column col for the sample"""
    table = get_table(wildcards.cohort)
    value = table.loc[wildcards.sample, col]
    return value


def get_tumor_type_mskcc_oncotree(wildcards):
    """Get the tumor type MSKCC oncotree of the sample"""
    return get_column_table_sample(wildcards, "MSKCC_Oncotree")


def get_tumor_type_civic(wildcards):
    """Get the tumor type Civic_Disease of the sample"""
    return get_column_table_sample(wildcards, "Civic_Disease")
