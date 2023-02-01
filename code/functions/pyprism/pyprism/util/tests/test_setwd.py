# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 2020

@author: Yoann Pradat

Tests for _setwd.py module.
"""

import os
from pyprism.util import setwd_to_data, setwd_to_results

def test_setwd_data():
    cwd = setwd_to_data()
    assert os.getcwd().endswith("data")

    dat_folders = ["cln", "rna", "wes", "util"]
    cwd_folders = os.listdir(os.getcwd())
    assert set(cwd_folders).issubset(set(cwd_folders))

    os.chdir(cwd)

def test_setwd_results():
    cwd = setwd_to_results()
    assert os.getcwd().endswith("results")
    os.chdir(cwd)
