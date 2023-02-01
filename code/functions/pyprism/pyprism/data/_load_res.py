# -*- coding: utf-8 -*-
"""
@created: Jul 08 2021
@modified: Jul 08 2021
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Function for load gene lists.
    load_gene_list
"""

import os
import pandas as pd

from ._util_load import load_from_data
from ._filepaths import _get_filepath_res

DataFrame = pd.core.frame.DataFrame

# User function ========================================================================================================

def load_resource(database, name, **kwargs):
    """
    Load lists of genes from different databases or used for different analyses.

    For all sources, the gene symbol is encoded in the variable \code{Hugo_Symbol}. Additional gene identifiers include
    the Ensembl gene id \code{Ensembl_Gene_Id} and the Entrez gene id \code{Entrez_Gene_Id}.

    Parameters
    ----------
    database: str
        Database representing the resource. Supported values are
        - civic
        - cosmic
        - gencode
        - oncokb
    name: str
        Name of the resource from the database. Supported values are
        - curated -> gene
        - cosmic -> gene
        - civic -> assertion, clinical_evidence, gene, variant, variant_group
        - oncokb -> gene
        - gencode -> gencode_v23, gencode_v27

    Returns
    -------
    df: DataFrame
        A dataframe with the gene list loaded and possible extra annotations of genes.
    """
    df = load_from_data(_get_filepath_res(database.lower(), name.lower()), **kwargs)

    return df
