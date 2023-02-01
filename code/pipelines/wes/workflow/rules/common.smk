from itertools import product
import pandas as pd
import re
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.4.0")

B_FOLDER = "workflow/benchmarks"
L_FOLDER = "workflow/logs"
R_FOLDER = "results"

###### Config file and sample sheets #####
configfile: "config/config.yaml"

table = pd.read_table(config["samples"], dtype=str).set_index(["Sample_Id"], drop=False)
samples = table["Sample_Id"].tolist()
nsamples = table.loc[(table["Sample_Type"]=="DNA_N")]["Sample_Id"].tolist()
nsamples_na = nsamples + ["NA"]
tsamples = table.loc[(table["Sample_Type"]=="DNA_T")]["Sample_Id"].tolist()
fastqs = ["%s_R%d" % (sample, i) for sample in samples for i in [1,2]]

if "normal_pool" in config and config["normal_pool"] is not None:
    table_normal_pool = pd.read_table(config["normal_pool"], dtype=str)
    nsamples_normal_pool = table_normal_pool["Sample_Id"].tolist()
    nsample_normal_pool = "_".join(table_normal_pool["Sample_Id"].tolist())
else:
    nsamples_normal_pool = []
    nsample_normal_pool = ""

def get_intervals_for_pon_target_file(target):
    folder = config["target_files"]["intervals_pon"][target]
    file = "%s_{interval}.bed" % config["target_files"]["target"][target]
    intervals, = glob_wildcards(os.path.join(folder, file))
    return intervals


if "Capture_Kit" in table:
    targets = table["Capture_Kit"].unique().tolist()
    intervals = set().union(*[get_intervals_for_pon_target_file(target) for target in targets])
else:
    targets = []
    intervals = []

if "download_gdc" in config.keys():
    ids_gdc_tcga_pancanatlas = pd.read_table(config["download_gdc"]["manifests"]["tcga_pancanatlas"])["id"].tolist()


##### Helper functions #####

def filter_combinator(combinator, comblist, white_list=True):
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
        # Use frozenset instead of tuple
        # in order to accomodate
        # unpredictable wildcard order
            if white_list:
                if frozenset(wc_comb) in comblist:
                    yield wc_comb
            else:
                if frozenset(wc_comb) not in comblist:
                    yield wc_comb
    return filtered_combinator


def get_allowed_pairs_tumor_normal():
    allowed = []
    df_pairs = pd.read_table(config["tumor_normal_pairs"]).fillna("NA")
    for (tsample, nsample) in zip(df_pairs["DNA_T"], df_pairs["DNA_N"]):
        allowed.append(frozenset({("tsample", tsample), ("nsample", nsample)}))
    return filter_combinator(product, allowed, white_list=True)


def get_allowed_pairs_target_interval():
    allowed = []
    for target in targets:
        folder = config["target_files"]["intervals_pon"][target]
        file = "%s_{interval}.bed" % config["target_files"]["target"][target]
        intervals, = glob_wildcards(os.path.join(folder, file))
        for interval in intervals:
            allowed.append(frozenset({("target", target), ("interval", interval)}))
    return filter_combinator(product, allowed, white_list=True)


def get_fastqs_irods(wildcards):
    """Get fastq files irods paths of given sample."""
    fastqs = table.loc[wildcards.sample, ["FASTQ_1", "FASTQ_2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.FASTQ_1, "r2": fastqs.FASTQ_2}
    return {"r1": fastqs.FASTQ_1}


def get_fastqs_local(wildcards):
    """Get fastq files local paths of given sample."""
    return {"r%d" % i: "%s/data/fastq/%s_R%d.fastq.gz" % (R_FOLDER, wildcards.sample, i) for i in [1,2]}


def get_fastqs_names(wildcards):
    """Get fastq files base name of given sample."""
    fastqs = table.loc[wildcards.sample, ["FASTQ_1_Name", "FASTQ_2_Name"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.FASTQ_1_Name, "r2": fastqs.FASTQ_2_Name}
    return {"r1": fastqs.FASTQ_1_Name}


def get_bam_id(wildcards):
    """Get bam id of given bam name."""
    return table.loc[table["BAM_Name"]==wildcards.bam_name]["BAM_Id"].iloc[0]


def get_column_table_sample(wildcards, col):
    """Get the value of the column col for the sample"""
    try:
        value = table.loc[wildcards.sample, col]
    except AttributeError:
        try:
            value = table.loc[wildcards.tsample, col]
        except AttributeError:
            try:
                value = table.loc[wildcards.nsample, col]
            except AttributeError:
                if wildcards.sample_pair=="all_samples":
                    value = ""
                else:
                    tsample = wildcards.sample_pair.split("_vs_")[0]
                    value = table.loc[tsample, col]
    return value


def get_target_name_sample(wildcards):
    """Get name of target file used for the sample"""
    return get_column_table_sample(wildcards, "Capture_Kit")


def get_target_file_sample(wildcards, file="bed_padded"):
    """Get path to target files used for the sample"""
    capture_kit = get_target_name_sample(wildcards)
    return config["target_files"][file][capture_kit]


def get_target_name_tumor_normal(wildcards):
    """Get name of target file used for the pair tumor, normal"""
    try:
        n_capture_kit = table.loc[wildcards.nsample, "Capture_Kit"]
        t_capture_kit = table.loc[wildcards.tsample, "Capture_Kit"]
        if n_capture_kit!=t_capture_kit:
            capture_kit = [n_capture_kit, t_capture_kit]
            intersect = True
        else:
            capture_kit = t_capture_kit
            intersect = False
    except AttributeError:
        t_capture_kit = table.loc[wildcards.tsample, "Capture_Kit"]
        capture_kit = t_capture_kit
        intersect = False
    return capture_kit, intersect


def get_target_file_tumor_normal(wildcards, file="bed_padded"):
    """Get path to target file to be used when calling somatic mutations.
    It may happen that the targets used for the normal and the tumor are different.
    In this case, the intersection of all target files of the project is used
    as the target file.
    """
    capture_kit, intersect = get_target_name_tumor_normal(wildcards)
    if intersect:
        folder = config["target_files"][file]["intersections"]
        capture_kit_0 = config["target_files"]["target"][capture_kit[0]]
        capture_kit_1 = config["target_files"]["target"][capture_kit[1]]

        if file=="bed":
            fn_1 = "%s_vs_%s.bed" % (capture_kit_0, capture_kit_1)
            fn_2 = "%s_vs_%s.bed" % (capture_kit_1, capture_kit_0)
            fp_1 = os.path.join(folder, fn_1)
            fp_2 = os.path.join(folder, fn_2)
            if os.path.exists(fp_1):
                target_file = fp_1
            elif os.path.exists(fp_2):
                target_file = fp_2
            else:
                raise ValueError("none of %s and %s exist!" % (fp_1, fp_2))

        elif file=="bed_padded":
            fn_1 = "%s_vs_%s_padded_10n.bed" % (capture_kit_0, capture_kit_1)
            fn_2 = "%s_vs_%s_padded_10n.bed" % (capture_kit_1, capture_kit_0)
            fp_1 = os.path.join(folder, fn_1)
            fp_2 = os.path.join(folder, fn_2)
            if os.path.exists(fp_1):
                target_file = fp_1
            elif os.path.exists(fp_2):
                target_file = fp_2
            else:
                raise ValueError("none of %s and %s exist!" % (fp_1, fp_2))

        else:
            raise ValueError("intersection files for %s have not yet been generated" % file)

        return target_file
    else:
        return config["target_files"][file][capture_kit]


def get_samples_target_file(wildcards, type="DNA_T"):
    """Return the list of samples prepared with a given target file."""
    samples_sel = []
    for sample, sample_type, sample_target in zip(table["Sample_Id"], table["Sample_Type"], table["Capture_Kit"]):
        if sample_type==type and sample_target==wildcards.target:
            samples_sel.append(sample)
    return samples_sel


def get_input_vcfs_pon_genomicsdb(wildcards):
    nsamples_target = get_samples_target_file(wildcards, type="DNA_N")
    return expand("%s/calling/pon_mutect2/{nsample}.vcf.gz" % R_FOLDER, nsample=nsamples_target)


def get_input_bed_pon_genomicsdb(wildcards):
    folder = config["target_files"]["intervals_pon"][wildcards.target]
    file = "%s_%s.bed" % (config["target_files"]["target"][wildcards.target], wildcards.interval)
    return os.path.join(folder, file)


def get_annovar_operation():
    protocols=config["params"]["annovar"]["protocol"]
    operations=[]
    for protocol in protocols:
        if protocol.startswith("refGene"):
            operations.append("g")
        else:
            operations.append("f")
    return ",".join(operations)


def get_target_names():
	return sorted(list(config["target_files"]["target"].keys()))


def get_target_files(names, file="bed"):
    return [config["target_files"][file][name] for name in names]


def get_tumor_type_mskcc_oncotree(wildcards):
    """Get the tumor type MSKCC oncotree of the sample"""
    return get_column_table_sample(wildcards, "MSKCC_Oncotree")


def get_tumor_type_civic(wildcards):
    """Get the tumor type Civic_Disease of the sample"""
    return get_column_table_sample(wildcards, "Civic_Disease")


def get_facets_diplogr(wildcards):
    """Get manually chosen value of dipLogR for FACETS, if any."""
    try:
        df_diplogr = pd.read_table(config["params"]["cnv"]["cnv_facets"]["diplogr"]).set_index(["DNA_P"])
        try:
            dna_p = "%s_vs_%s" % (wildcards.tsample, wildcards.nsample)
        except:
            dna_p = "%s_vs_NA" % wildcards.tsample
        diplogr = df_diplogr.loc[dna_p, "DipLogR"]
        return diplogr
    except:
        return -999


def get_input_somatic_maf_filter_false_positives(w):
    path_1 = "%s/calling/somatic_maf_filter_mutect_calls/%s_vs_%s.vcf.gz" % (R_FOLDER, w.tsample, w.nsample)
    path_2 = "%s/calling/somatic_maf2vcf/%s_vs_%s.vcf.gz" % (R_FOLDER, w.tsample, w.nsample)
    if w.tsample.startswith("TCGA"):
        return path_2
    else:
        return path_1


def get_input_concatenate(w, typ, db):
    input_folder = "%s/annotation/somatic_%s_%s" % (R_FOLDER, typ, db)

    if config["params"][db]["run_per_sample"][typ]:
        sample_pairs = expand("{tsample}_vs_{nsample}", get_allowed_pairs_tumor_normal(),
            tsample=tsamples, nsample=nsamples_na)
    else:
        sample_pairs = ["all_samples"]

    if typ=="maf":
        return ["%s/%s.maf" % (input_folder, sample_pair) for sample_pair in sample_pairs]
    elif typ=="cna":
        return ["%s/%s.tsv" % (input_folder, sample_pair) for sample_pair in sample_pairs]
