from collections import OrderedDict
from itertools import product
import ntpath
import pandas as pd
from snakemake.utils import min_version
import re
import yaml

min_version("5.4.0")

def unnest_list(v):
    if type(v)==list and any([type(e) in [list, dict, OrderedDict] for e in v]):
        unnest = True
    else:
        unnest = False

    if not unnest:
        if type(v)!=list:
            if type(v)==dict or type(v)==OrderedDict:
                unnested = []
                for e in list(v.values()):
                    unnested += unnest_list(e)
                return unnested
            else:
                return [v]
        else:
            return v
    else:
        unnested = []
        for e in v:
            unnested += unnest_list(e)
        return unnested

configfile: "config/config.yaml"

L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
R_FOLDER = "../../../results/survival_analysis"
D_FOLDER = config["filepaths"]["base"]
F_FOLDER = "../../../results/figures_paper"
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)

COHORTS = list(config["data"].keys())
SAMPLES = list(config["models"]["samples"].keys())
FEATURES = list(config["models"]["features"].keys())
SELECTIONS = config["models"]["selections"]
MODELS = config["models"]["models"]

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


def get_allowed_doublets_samples_features():
    allowed = []
    for cohort in COHORTS:
        for samples in SAMPLES:
            for features in config["models"]["samples"][samples]:
                allowed.append(frozenset({("cohort", cohort), ("samples", samples), ("features", features)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_tuples_samples_features_selection_model():
    allowed = []
    config_features_models = config["models"]["features_models"]
    config_features_selections = config["models"]["features_selections"]

    for cohort in COHORTS:
        for samples in SAMPLES:
            for features in config["models"]["samples"][samples]:
                if config_features_models is not None and features in config_features_models:
                    models_all = config_features_models[features]
                else:
                    models_all =  MODELS

                for model in models_all:
                    if config_features_selections is not None and features in config_features_selections:
                        if config_features_selections[features] is not None:
                            selections_dict = [s for s in config_features_selections[features] if type(s)==dict \
                                               or type(s)==OrderedDict]
                            if any(model in d for d in selections_dict):
                                for d in selections_dict:
                                    if model in d:
                                        selections_all = d[model]
                                        break
                            else:
                                selections_all = config_features_selections[features]
                                selections_all = [s for s in selections_all if type(s)!=dict and type(s)!=OrderedDict]
                        else:
                            selections_all = config_features_selections[features]
                            selections_all = [s for s in selections_all if type(s)!=dict and type(s)!=OrderedDict]
                    else:
                        selections_all = SELECTIONS

                    for selection in selections_all:
                        allowed.append(frozenset({("cohort", cohort), ("samples", samples), ("features", features),
                                                  ("selection", selection), ("model", model)}))
    return filter_combinator(product, allowed, white_list=True)


def get_input_indices_run_survival(w):
    repeat_mode = config["models"]["fitting"]["repeat_mode"]
    return "%s/models_{cohort}/{samples}/{features}/%s_indices.tsv.gz" % (R_FOLDER, repeat_mode)


def get_input_pool_ax_imputations(w):
    checkpoint_output = checkpoints.group_encode_check_impute.get(**w).output[0]
    tables_checkpoint = glob_wildcards(os.path.join(checkpoint_output, "data.{table}.tsv.gz")).table
    folder = "%s/models_%s/%s/%s/%s_%s" % (R_FOLDER, w.cohort, w.samples, w.features, w.selection, w.model)
    return ["%s/%s.%s.tsv.gz" % (folder, data, table) for data in ["mets", "covs"] for table in tables_checkpoint]


def get_input_plot_metrics(w):
    config_features_models = config["models"]["features_models"]
    config_features_selections = config["models"]["features_selections"]
    features_all = config["models"]["samples"][w.samples]
    inputs = []

    for features in features_all:
        if config_features_models is not None and features in config_features_models:
            models_all = config_features_models[features]
        else:
            models_all =  MODELS

        for model in models_all:
            if config_features_selections is not None and features in config_features_selections:
                if config_features_selections[features] is not None:
                    selections_dict = [s for s in config_features_selections[features] if type(s)==dict \
                                       or type(s)==OrderedDict]
                    if any(model in d for d in selections_dict):
                        for d in selections_dict:
                            if model in d:
                                selections_all = d[model]
                                break
                    else:
                        selections_all = config_features_selections[features]
                        selections_all = [s for s in selections_all if type(s)!=dict and type(s)!=OrderedDict]
                else:
                    selections_all = config_features_selections[features]
                    selections_all = [s for s in selections_all if type(s)!=dict and type(s)!=OrderedDict]
            else:
                selections_all = SELECTIONS

            if w.selection in selections_all:
                folder = "%s/models_%s/%s/%s/%s_%s" % (R_FOLDER, w.cohort, w.samples, features, w.selection, model)
                inputs.append("%s/mets.pooled_ax_rep.tsv.gz" % folder)
                inputs.append("%s/covs.pooled_ax_rep.tsv.gz" % folder)

    return inputs


def get_input_dats_plot_discretized_risk_score(w):
    checkpoint_output = checkpoints.group_encode_check_impute.get(**w).output[0]
    tables_checkpoint = glob_wildcards(os.path.join(checkpoint_output, "data.{table}.tsv.gz")).table
    folder = "%s/data_%s/sub_features/%s/%s/processed" % (R_FOLDER, w.cohort, w.samples, w.features)
    return ["%s/data.%s.tsv.gz" % (folder, table) for table in tables_checkpoint]


def get_name_model_plot_discretized_risk_score(w):
    features2model = {"cln_grim_0": "M1", "cln_biol_1": "M2", "cln_biol_1_dna_1": "M3", "cln_biol_1_dna_2": "M4",
                      "cln_biol_1_rna_1": "M5", "cln_biol_1_rna_2": "M6", "cln_biol_1_dna_1_rna_1": "M7",
                      "cln_dna_1_rna_1": "M7bis"}

    if w.features in features2model:
        return features2model[w.features]
    else:
        return "MX"
