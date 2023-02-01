import numpy as np
import pandas as pd

def filter_regex(df, regex, column, mode="out", tag_only=False, tag=None):
    mask_regex = df[column].str.contains(regex)
    mask_regex = mask_regex.replace(np.nan, False)
    if mode=="out":
        print("-filtered out %d/%d fusions using column %s and regex %s" % \
              (sum(mask_regex), mask_regex.shape[0], column, regex))
    elif mode=="in":
        mask_regex = ~mask_regex
        print("-filtered in %d/%d fusions using column %s and regex %s" % \
              (sum(~mask_regex), mask_regex.shape[0], column, regex))
    else:
        raise ValueError("-ERROR: unrecognized value of mode %s" % mode)

    if tag_only:
        if "FILTER" not in df.columns:
            df.loc[mask_regex, "FILTER"] = tag
        else:
            mask_null = df["FILTER"].isnull()
            df.loc[(mask_regex) & (~mask_null), "FILTER"] += ";%s" % tag
            df.loc[(mask_regex) & (mask_null), "FILTER"] = tag
        return df
    else:
        return df.loc[~mask_regex]


def filter_out_potential_fp(df, column, mode="FusionAnnotator", tag_only=False):
    if mode=="FusionAnnotator":
        # list taken from the "Red Herrings" section of https://github.com/FusionAnnotator/CTAT_HumanFusionLib/wiki
        potential_fp = ["GTEx_recurrent_StarF2019", "BodyMap", "DGD_PARALOGS", "HGNC_GENEFAM", "Greger_Normal",
                        "Babiceanu_Normal", "ConjoinG"]
        tag = "Red_Herrings_FusionAnnotator"
    elif mode=="Custom":
        potential_fp = ["Babiceanu_Normal", "ChimerSeq_Normal_v4.0", "GTEX_V6_NAR_2020"]
        tag = "Found_In_Normal_Tissues"
    else:
        raise ValueError("-ERROR: unrecognized value of mode %s" % mode)

    regex_potential_fp = "|".join(potential_fp)
    return filter_regex(df, regex_potential_fp, column, mode="out", tag_only=tag_only, tag=tag)


def filter_out_gene_types(df, tag_only=False):
    # excluded gene types are taken from methods of PMCID: PMC5916809
    gene_types_fp = ["^TEC", "^rRNA", "^Mt_", "^IG_", "^TR_", "^noncoding", "^ribozyme", "^sRNA", "^scRNA", "^scaRNA",
                     "^pseudogene", "^NA"]
    regex_gene_types_fp = "|".join(gene_types_fp)
    df["Gene_Type_1"] = df["Gene_Type_1"].fillna("NA")
    df["Gene_Type_2"] = df["Gene_Type_2"].fillna("NA")
    tag_1 = "Excluded_Gene_1_Type"
    tag_2 = "Excluded_Gene_2_Type"
    df = filter_regex(df, regex_gene_types_fp, "Gene_Type_1", mode="out", tag_only=tag_only, tag=tag_1)
    df = filter_regex(df, regex_gene_types_fp, "Gene_Type_2", mode="out", tag_only=tag_only, tag=tag_2)

    return df


def filter_out_fusions(df, tag_only=False):
    df = filter_out_potential_fp(df, column="Annotations", mode="FusionAnnotator", tag_only=tag_only)
    df = filter_out_potential_fp(df, column="Annotations_Custom", mode="Custom", tag_only=tag_only)
    df = filter_out_gene_types(df, tag_only=tag_only)
    return df


def filter_in_whitelist(df, tag_only=False):
    column = "Annotations_Custom"
    whitelist = ["ChimerKB_v4.0", "Chitars_Cancer_v5.0", "COSMIC_curated_v95", "ONE_PARTNER_IS_DRIVER", "TIC_v3.3"]
    regex_whitelist = "|".join(whitelist)
    tag = "Not_In_Whitelist"
    return filter_regex(df, regex_whitelist, column, mode="in", tag_only=tag_only, tag=tag)
