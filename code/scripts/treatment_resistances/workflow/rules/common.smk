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
R_FOLDER = "../../../results/treatment_resistances"
F_FOLDER = "../../../results/figures_paper"
D_FOLDER = config["filepaths"]["base"]
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)

SAMPLES = list(config["models"]["samples"].keys())
FEATURES = list(config["models"]["features"].keys())
MODELS = config["models"]["models"]
OUTCOME_RGXS = list(set(unnest_list(list(config["models"]["outcome_rgxs"].values()))))
SELECTIONS = list(set(unnest_list(v=list(config["models"]["selections"].values()))))

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


def get_allowed_tuples_samples_features():
    allowed = []
    for samples in SAMPLES:
        for features in config["models"]["samples"][samples]:
            allowed.append(frozenset({("samples", samples), ("features", features)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_tuples_samples_features_outcome_rgx():
    allowed = []
    for samples in SAMPLES:

        samples_root = []
        for k,v in config["models"]["outcome_rgxs"].items():
            if k in samples:
                samples_root.append(k)
        if len(samples_root)!=1:
            raise ValueError("-0 or multiple matched sample root")
        samples_root = samples_root[0]

        for features in config["models"]["samples"][samples]:
            for outcome_rgx in config["models"]["outcome_rgxs"][samples_root]:
                allowed.append(frozenset({("samples", samples), ("features", features), ("outcome_rgx", outcome_rgx)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_tuples_samples_features_outcome_rgx_selection_model():
    allowed = []
    for samples in SAMPLES:
        samples_root = []
        for k,v in config["models"]["outcome_rgxs"].items():
            if k in samples:
                samples_root.append(k)
        if len(samples_root)!=1:
            raise ValueError("-0 or multiple matched sample root")
        samples_root = samples_root[0]

        for features in config["models"]["samples"][samples]:
            if features in config["models"]["selections"]:
                selections_features = config["models"]["selections"][features]
            else:
                selections_features = config["models"]["selections"]["default"]

            for outcome_rgx in config["models"]["outcome_rgxs"][samples_root]:
                for model in MODELS:
                    if model in selections_features:
                        selections = selections_features[model]
                    else:
                        selections_all = selections_features
                        selections = [s for s in selections_all if type(s)!=dict and type(s)!=OrderedDict]

                    for selection in selections:
                        allowed.append(frozenset({("samples", samples), ("features", features),
                                                  ("outcome_rgx", outcome_rgx), ("model", model),
                                                  ("selection", selection)}))
    return filter_combinator(product, allowed, white_list=True)


def get_input_pool_ax_imputations(w):
    checkpoint_output = checkpoints.group_encode_check_impute.get(**w).output[0]
    tables_checkpoint = glob_wildcards(os.path.join(checkpoint_output, "data.{table}.tsv.gz")).table
    folder = "%s/models/%s/%s/%s/%s_%s" % (R_FOLDER, w.samples, w.features, w.outcome_rgx,  w.selection, w.model)
    return ["%s/%s.%s.tsv.gz" % (folder, data, table) for data in ["mets", "covs"] for table in tables_checkpoint]
