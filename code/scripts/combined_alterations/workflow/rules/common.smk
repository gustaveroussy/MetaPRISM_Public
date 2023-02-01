from snakemake.utils import min_version
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

R_FOLDER = "../../../results/combined_alterations"
F_FOLDER = "../../../results/figures_paper"
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
D_FOLDER = config["filepaths"]["base"]
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)
