from itertools import product
import ntpath
import pandas as pd
from snakemake.utils import min_version
import re
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
R_FOLDER = "../../../results/data_overview"
F_FOLDER = "../../../results/figures_paper"
D_FOLDER = config["filepaths"]["base"]
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)
