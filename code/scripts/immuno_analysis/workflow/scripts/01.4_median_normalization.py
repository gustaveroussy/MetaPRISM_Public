# -*- coding: utf-8 -*-
"""
Created: 28/07/21
Modified: 28/07/21
Author: Antoine Lainé
"""

import argparse
import sys
import csv
import pandas as pd
import os

## Functions ##
def normalize_signature(RawSignature,Annotation):
    #Split by Project_TCGA / Divide each part by median / Join new tables
    Annotation.sort_values(by="Project_TCGA", axis=0, inplace=True)
    Annotation.set_index(keys=['Project_TCGA'], drop=False,inplace=True)
    names=Annotation['Project_TCGA'].unique().tolist()
    NormalizedSignaData = []
    for i in names:
        row_iter = Annotation.loc[Annotation.Project_TCGA==i]
        samples = row_iter.Sample_Id
        sub = RawSignature[samples]
        normalized = sub.div(sub.median(axis=1), axis=0)
        try:
            NormalizedSignaData = pd.concat([NormalizedSignaData, normalized], axis=1)
        except:
            NormalizedSignaData = normalized
    NormalizedSignature = NormalizedSignaData
    return NormalizedSignature


## Parameters ##

parser = argparse.ArgumentParser()
parser.add_argument('--output', type=str ,help='Output folder')
parser.add_argument('--selected_path', type=str,help='Path to selected_samples.tsv')
parser.add_argument('--raw_signature_path',type=str ,help='Path to raw signature')
args = parser.parse_args()

raw_signature_bn = os.path.basename(args.raw_signature_path).replace(".tsv", "")
## Main ##

raw_signature = pd.read_csv(args.raw_signature_path,sep="\t",index_col=0)
annotation = pd.read_csv(args.selected_path,sep="\t")
normalized_signature = normalize_signature(raw_signature,annotation)
normalized_signature.to_csv(args.output + "/" + raw_signature_bn + "_normalized.tsv",sep="\t")
