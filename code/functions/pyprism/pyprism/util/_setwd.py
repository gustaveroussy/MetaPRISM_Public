# -*- coding: utf-8 -*-
"""
@modified: 06 Nov 2020
@created: 30 Sep 2020
@author: Yoann Pradat

    CentraleSupelec
    MICS laboratory
    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France

    Institut Gustave Roussy
    Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Fonctions for setting the current working directory to folders of the project.

The following functions are implemented
    setwd_to_data()
    setwd_to_results()
    setwd_to_scripts()
"""

import os

# Name of the project root folder. Modify if you've cloned the git repository under a different name
PROJECT_FOLDER = "MetaPRISM"

def setwd_to_data():
    current_wd = os.getcwd()
    if PROJECT_FOLDER not in os.getcwd():
        raise ValueError("Please set the working directory to a location in the repository %s" % PROJECT_FOLDER)
    else:
        while not os.getcwd().endswith(PROJECT_FOLDER):
            os.chdir("..")
        os.chdir("./data/")
    return current_wd

def setwd_to_results():
    current_wd = os.getcwd()
    if PROJECT_FOLDER not in os.getcwd():
        raise ValueError("Please set the working directory to a location in the repository %s" % PROJECT_FOLDER)
    else:
        while not os.getcwd().endswith(PROJECT_FOLDER):
            os.chdir("..")
        os.chdir("./results/")
    return current_wd


def setwd_to_scripts():
    current_wd = os.getcwd()
    if PROJECT_FOLDER not in os.getcwd():
        raise ValueError("please set the working directory to a location in the repository %s" % PROJECT_FOLDER)
    else:
        while not os.getcwd().endswith(PROJECT_FOLDER):
            os.chdir("..")
        os.chdir("./scripts/")
    return current_wd
