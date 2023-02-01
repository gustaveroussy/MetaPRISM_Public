import ntpath
import pandas as pd
from snakemake.utils import min_version
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

D_FOLDER = "../../../data"
R_FOLDER = "../../../results/mutational_signatures"
F_FOLDER = "../../../results/figures_paper"
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)
