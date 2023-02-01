from snakemake.utils import min_version
import yaml

min_version("5.4.0")

configfile: "config/config.yaml"

R_FOLDER = "../../../results/germline_mutations"
F_FOLDER = "../../../results/figures_paper"
L_FOLDER = "workflow/logs"
B_FOLDER = "workflow/benchmarks"
D_FOLDER = config["filepaths"]["base"]
FILEPATHS = yaml.load(open(config["filepaths"]["yaml"], "r"), Loader=yaml.FullLoader)

#### Helper functions ####
def get_mutations_file(w):
    if w.selection_mut == "annotated":
        return config["data"]["mutations"]["ann"][w.cohort]
    else:
        return config["data"]["mutations"]["all"][w.cohort]
