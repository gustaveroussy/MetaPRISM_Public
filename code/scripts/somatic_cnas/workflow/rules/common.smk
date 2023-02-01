from snakemake.utils import min_version
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

D_FOLDER = "../../../data"
F_FOLDER = "../../../results/figures_paper"
R_FOLDER = "../../../results/somatic_cnas"
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)
