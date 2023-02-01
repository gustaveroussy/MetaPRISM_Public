import ntpath
import pandas as pd
from snakemake.utils import min_version
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

R_FOLDER = "../../../results/immuno_analysis"
F_FOLDER = "../../../results/figures_paper"
D_FOLDER = config["filepaths"]["base"]
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)

DATA_COHORTS = config["data"]["cohorts"]

MFP_MODELS_COHORTS = config["mfp_models"]["cohorts"]
MFP_MODELS_NAMES = config["mfp_models"]["names"]
MFP_PREDICT_COHORTS = list(set(DATA_COHORTS).union(set(MFP_MODELS_COHORTS)))

##### Helper functions #####

def get_gene_table_udpated(w):
    gene_table = config["resources"][w.cohort]["gene_table"]
    return  "%s/gene_tables/%s_updated.tsv" % (R_FOLDER, gene_table)


def get_gene_selection_udpated(w):
    gene_selection = config["resources"][w.cohort]["gene_selection"]
    return  "%s/gene_tables/%s_updated.tsv" % (R_FOLDER, gene_selection)


def get_expression_data(w):
    path_1 = D_FOLDER
    path_2 = FILEPATHS[w.cohort]["rna_gex"]["full"]["genes"]["tpm"]["data"]
    return os.path.join(path_1, path_2)


def get_expression_summary(w):
    path_1 = D_FOLDER
    path_2 = FILEPATHS[w.cohort]["rna_gex"]["full"]["genes"]["tpm"]["summary"]
    return os.path.join(path_1, path_2)

# MFP model training

def get_mfp_model_signatures(w):
    if w.cohort == "bagaev_2021":
        return config["resources"]["bagaev_2021"]["signatures"]
    elif w.cohort == "tcga":
        return "%s/%s/tables/signatures_mfp_model.tsv" % (R_FOLDER, w.cohort)


def get_mfp_model_annotation(w):
    if w.cohort == "bagaev_2021":
        return config["resources"]["bagaev_2021"]["annotation"]
    else:
        return "%s/%s/tables/annotation_mfp_model.tsv" % (R_FOLDER, w.cohort)


def get_mem_mb_mfp_model(w):
    return config["mfp_models"][w.model_name]["mem_mb_per_thread"]*config["mfp_models"][w.model_name]["threads"]


def get_nodes_mfp_model(w):
    return config["mfp_models"][w.model_name]["nodes"]


def get_partition_mfp_model(w):
    return config["mfp_models"][w.model_name]["partition"]


def get_threads_mfp_model(w):
    return config["mfp_models"][w.model_name]["threads"]


def get_time_mfp_model(w):
    return config["mfp_models"][w.model_name]["time"]


def get_input_models_mfp_model(w):
    models_meta = [x for x in MFP_MODELS_NAMES if x.startswith("MetaA") or x.startswith("MetaB")]
    models_indv = [x for x in MFP_MODELS_NAMES if x not in models_meta]
    if w.model_name in models_indv:
        return []
    else:
        models = config["mfp_models"][w.model_name]["estimators"]
        return ["%s/%s/mfp_models/%s/pipeline_train.joblib" % (R_FOLDER, w.cohort, model) for model in models]

# MFP model prediction

def get_mfp_predict_samples(w):
    if w.cohort == "bagaev_2021":
        return []
    else:
        return "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER


def get_mfp_predict_signatures(w):
    if w.cohort == "bagaev_2021":
        return config["resources"]["bagaev_2021"]["signatures"]
    else:
        return "%s/%s/tables/signatures_bagaev_2021_mfp.tsv" % (R_FOLDER, w.cohort)


def get_mfp_predict_annotation(w):
    if w.cohort == "bagaev_2021":
        return config["resources"]["bagaev_2021"]["annotation"]
    else:
        return "pyprism_%s" % w.cohort


def get_mfp_predict_model(w):
    if w.cohort == "bagaev_2021":
        return "%s/bagaev_2021/mfp_models/%s/pipeline_train.joblib" % (R_FOLDER, w.model_name)
    else:
        return "%s/tcga/mfp_models/%s/pipeline_train.joblib" % (R_FOLDER, w.model_name)


def get_mfp_extract_model(w):
    if w.cohort == "bagaev_2021":
        return "%s/bagaev_2021/mfp_models/%s/pipeline_train.joblib" % (R_FOLDER, MFP_MODELS_NAMES[0])
    else:
        return "%s/tcga/mfp_models/%s/pipeline_train.joblib" % (R_FOLDER, MFP_MODELS_NAMES[0])
