import pandas

from collections import defaultdict
from sys import version_info
from snakemake.utils import min_version
from snakemake.remote.HTTP import RemoteProvider
from typing import Any, Dict

HTTP = RemoteProvider()

#########################
### ENVS REQUIREMENTS ###
#########################
if version_info < (3, 8):
    raise SystemError("Please use Python 3.8 or later")

min_version('5.20.1')
container: "docker://continuumio/miniconda3:4.4.10"
git = "https://raw.githubusercontent.com/tdayris/snakemake-wrappers/892df8bb84ec95f027ab5ff0a65787902cb65fd1"

if config is None or config == {}:
    config: "config.yaml"

design = pandas.read_csv(
    config["design"],
    sep="\t",
    header=0,
    index_col=None,
    dtype=str
)

if any(design.Sample_id.duplicated()):
    raise ValueError("Sample ID must be unique")

########################
### GLOBAL VARIABLES ###
########################

def get_fasta(config: Dict[str, Any]):
    """
    According to parameters, return the correct fasta file
    """
    if config["database"].get("db_name", "gencode") == "gencode":
        if config["database"].get("content", "all").lower() == "all":
            return "references/gencode.v27.transcripts.fa"
        elif config["database"].get("content", "all").lower() == "lncrna":
            return "resources/gencode.v27.lncRNA.transcripts.fa"
        elif config["database"].get("content", "all").lower() == "proteincoding":
            return "resources/gencode.v27.proteincoding.transcripts.fa"
        else:
            raise ValueError("Could not determine fasta file")
    elif config["database"].get("db_name", "gencode") == "NONCODE":
        return "resources/GENCODE_NONCODE.transcripts.fa"
    else:
        raise ValueError("Could not determine fasta file")


fasta_path = get_fasta(config)
gtf_path = (
    "references/gencode.v27.annotation.gtf"
    if config["database"].get("db_name", "gencode")
    else "resources/GENCODE_NONCODE.annotation.gtf"
)

############
### URLs ###
############
# This corresponds to the lines:
# https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh
# 13 to 15, 36 to 39
ref = {
    "gencode.v27.transcripts.fa": "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz",
    "gencode.v27.annotation.gtf": "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz",
    "GRCh38.primary_assembly.genome.fa": "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz",
    "NONCODEv5_human_hg38_lncRNA.gtf": "http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz",
    "NONCODEv5_human.fa": "http://www.noncode.org/datadownload/NONCODEv5_human.fa.gz"
}


# This links the scripts and their URL
scripts = {
    "Extract.pl": "https://raw.githubusercontent.com/gevaertlab/RNASeq_pipeline/master/scripts/Extract.pl",
    "Exclude.pl": "https://raw.githubusercontent.com/gevaertlab/RNASeq_pipeline/master/scripts/Exclude.pl",
    "match.pl": "https://raw.githubusercontent.com/gevaertlab/RNASeq_pipeline/master/scripts/match.pl"
}

# Additional resources provided
resources = {
    "duplicate_clusters.tsv": "https://raw.githubusercontent.com/gevaertlab/RNASeq_pipeline/master/prepare_reference/duplicate_clusters.tsv"
}

####################
### FILE PAIRING ###
####################

def get_subread_pairs(wildcards: Any) -> Dict[str, str]:
    """
    Return subread pairs
    """
    pairs = {
        "NONCODEv5_human_hg38_lncRNA": {
            "gtf": "resources/GENCODE_NONCODE.annotation.gtf",
            "fasta": "resources/GENCODE_NONCODE.transcripts.fa"
        },
        "gencode.v27": {
            "gtf": "references/gencode.v27.annotation.gtf",
            "fasta": "references/GRCh38.primary_assembly.genome.fa"
        }
    }

    return pairs[wildcards.database]


copy_dict = defaultdict(None)
for file, sample in zip(design.Upstream_file, design.Sample_id):
    copy_dict[f"{sample}_1.fastq.gz"] = file
for file, sample in zip(design.Downstream_file, design.Sample_id):
    copy_dict[f"{sample}_2.fastq.gz"] = file


def is_paired(wildcards) -> bool:
    """
    Return wether a sample is pair or single ended
    """
    try:
        return os.path.exists(copy_dict[f"{wildcards.sample}_2.fastq.gz"])
    except TypeError:
        return False
    #return os.path.exists(f"raw_data/{wildcards.sample}_2.fastq.gz")


def get_input_fq_files_kallisto(wildcards):
    """
    Return input files for kallisto. If sample is pair-ended, then it will be a
    list of two files, else the list will contain only one file
    """
    if is_paired(wildcards):
        return {"fastq": [
            f"trim_galore/fastq/{wildcards.sample}_1_val_1.fq.gz",
            f"trim_galore/fastq/{wildcards.sample}_2_val_2.fq.gz"
        ]}
    return {
        "fastq": [f"trim_galore/fastq/{wildcards.sample}_1_trimmed.fq.gz"],
        "fragment_size": f"fragment_size/{wildcards.sample}.txt"
    }
##################################
### SNAKEMAKE'S TARGET CALLING ###
##################################

wildcard_constraints:
    sample = "|".join(design.Sample_id),
    gencode = "|".join(ref.keys()),
    subset = "|".join(["proteincoding", "lncRNA"]),
    script = "|".join(scripts.keys())

localrules:
    copy_fastq


rule all:
    input:
        tximport = "tximport/txi.rds",
        kallisto = expand(
            "kallisto/quant/{sample}.tar.bz2",
            sample=design["Sample_id"].tolist()
        ),
        quantification = "tximport/quantification.csv"
    message:
        "Finishing pipeline"


#####################
### CORE PIPELINE ###
#####################

#################
### REFERENCE ###
#################
"""
This rule corresponds to the bash line:
https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh

13 to 15, 36 and 38
"""
rule download_reference:
    input:
        lambda wildcards: HTTP.remote(ref[wildcards.gencode])
    output:
        temp("references/{gencode}.gz")
    message:
        "Downloading {wildcards.gencode}.gz reference from ebi"
    group:
        "Download_Reference"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 2048),
        time_min = lambda wildcards, attempt: attempt * 120
    conda:
        "../../envs/bash_perl.yaml"
    params:
        extra = "--verbose"
    log:
        "logs/download/{gencode}.log"
    shell:
        "mv {params.extra} {input} {output} 2> {log}"


"""
This rule corresponds to the bash line:
https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh

Lines 16, 37, 39
"""
rule unzip_reference:
    input:
        "references/{gencode}.gz"
    output:
        temp("references/{gencode}")
    group:
        "Download_Reference"
    message:
        "Unzipping {wildcards.gencode}.gz"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 2048),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    log:
        "logs/gunzip/{gencode}.log"
    shell:
        "gunzip -c {input} > {output} 2> {log}"


"""
This rule corresponds to the bash line:
https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh#L20
"""
rule gene_transcript_type:
    input:
        "references/gencode.v27.annotation.gtf"
    output:
        temp("subsets/gene_transcript_type")
    message:
        "Building gene transcript type"
    threads:
        8
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 2048),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "envs.bash.yaml"
    log:
        less0 = "logs/gene_transcript_type/less0.log",
        grep1 = "logs/gene_transcript_type/grep1.log",
        awk2 = "logs/gene_transcript_type/awk2.log",
        cut3 = "logs/gene_transcript_type/cut3.log",
        awk4 = "logs/gene_transcript_type/awk4.log",
        sed5 = "logs/gene_transcript_type/sed5.log",
        awk6 = "logs/gene_transcript_type/awk6.log",
        sed7 = "logs/gene_transcript_type/sed7.log"
    shell:
        """less {input} 2> {log.less0} | grep -v ^# 2> {log.grep1} | awk '$3~/transcript/' 2> {log.awk2} | cut -f9 2> {log.cut3} | awk 'OFS=\"\\t\"{{print $2,$4,$6,$12}}' 2> {log.awk4} | sed 's/[";]//g' 2> {log.sed5} | awk 'OFS="\t"{{$5="other";$6="other";if($4~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$6="lncRNA";if($4~/protein_coding/)$6="proteincoding";if($3~/protein_coding/)$5="proteincoding";if($3~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$5="lncRNA";  print $0}}' 2> {log.awk6} | sed '1i gene_id\ttranscript\tgenetype1\ttranscripttype1\tgenetype\ttranscripttype' > {output} 2> {log.sed7}"""

"""
This rule corresponds to the bash line:
https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh#L21
"""
rule gene_type:
    input:
        "subsets/gene_transcript_type"
    output:
        temp("subsets/gene_type")
    message:
        "Building gene_type file"
    threads:
        6
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 2048),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "envs.bash.yaml"
    log:
        less0 = "logs/gene_type/less0.log",
        sed1 = "logs/gene_type/sed1.log",
        cut2 = "logs/gene_type/cut2.log",
        sort3 = "logs/gene_type/sort3.log",
        uniq4 =  "logs/gene_type/uniq4.log",
        sed5 =  "logs/gene_type/sed5.log"
    shell:
        "less {input} 2> {log.less0} | sed 1d 2> {log.sed1} | cut -f1,3,5 2> {log.cut2} | sort 2> {log.sort3} | uniq 2> {log.uniq4} | sed '1i gene_id\tdetail\ttype' > {output} 2> {log.sed5}"


"""
This rule corresponds to the bash lines:
https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh

22 and 23
"""
rule specific_gene_transcript_type:
    input:
        "subsets/gene_transcript_type"
    output:
        "subsets/gene_transcript_type.{subset}"
    message:
        "Building {wildcards.subset} specific gene_transcript_type"
    threads:
        2
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 2048),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/specific_gene_transcript_type/{subset}/less0.log",
        awk1 = "logs/specific_gene_transcript_type/{subset}/akw1.log",
    shell:
        "less {input} 2> {log.less0} | "
        "awk '$5==\"{wildcards.subset}\"' > {output} 2> {log.awk1}"


"""
This rule downloads all required scripts
"""
rule download_scripts:
    input:
        lambda wildcards: HTTP.remote(scripts[wildcards.script])
    output:
        "scripts/{script}"
    message:
        "Downloading {wildcards.script} from github"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 2048),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        "logs/download/{script}.log"
    shell:
        "wget {input} 2> {log}"


"""
This rule covers:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh

Line 25
"""
rule format_transcripts:
    input:
        "resources/gencode.v27.transcripts.fa"
    output:
        temp("resources/gencode.v27.transcripts.fa.tmp")
    message:
        "Format transcript fasta for further extractions"
    threads:
        5
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 2048),
        time_min = lambda wildcards, attempt: attempt * 45
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/format_transcripts/less0.log",
        awk1 = "logs/format_transcripts/awk1.log",
        sed2 = "logs/format_transcripts/sed2.log",
        sed3 = "logs/format_transcripts/sed3.log",
        sed4 = "logs/format_transcripts/sed4.log"
    shell:
        """less {input} 2> {log.less0} | awk '!/^>/ {{ printf "%s", $0; n = "\n" }} /^>/ {{ print n $0; n = "" }} END {{ printf "%s", n }}' 2> {log.awk1} | sed 's/^/%/' 2> {log.sed2} | sed '$!N;s/\n/\t/g' 2> {log.sed3} | sed 's/^%//' > {output} 2> {log.sed4}"""


"""
This rule covers:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
Lines 26 and 27
"""
rule surbset_transcripts_reference:
    input:
        fasta_tmp = "resources/gencode.v27.transcripts.fa.tmp",
        script = "scripts/Extract.pl",
        subset = "subsets/gene_transcript_type.{subset}"
    output:
        temp("subsets/gencode.v27.{subset}.transcripts.fa")
    message:
        "Subsetting fasta transcript on {wildcards.subset}"
    threads:
        6
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 8192),
        time_min = lambda wildcards, attempt: attempt * 45
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/surbset_transcripts_reference/{subset}/less0.log",
        awk1 = "logs/surbset_transcripts_reference/{subset}/awk1.log",
        sed2 = "logs/surbset_transcripts_reference/{subset}/sed2.log",
        extract3 = "logs/surbset_transcripts_reference/{subset}/extract3.log",
        cut4 = "logs/surbset_transcripts_reference/{subset}/cut4.log",
        sed5 = "logs/surbset_transcripts_reference/{subset}/sed5.log"
    shell:
        """less {input.fasta_tmp} 2> {log.less0} | awk -F "|" 'OFS="\t"{{print $1,$0}}' 2> {log.awk1} | sed 's/>//' 2> {log.sed2} | {input.script} {input.subset} 2 - 1 2> {log.extract3} | cut -f2- 2> {log.cut4} | sed 's/%/\n/'  > {output} 2> {log.sed5}"""


"""
This rule covers lines 29 and 30 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule subset_gtf_reference:
    input:
        gtf = "resources/gencode.v27.annotation.gtf",
        script = "scripts/Extract.pl",
        subset = "subsets/gene_transcript_type.{subset}"
    output:
        temp("subset/gencode.v27.{subset}.annotation.gtf")
    message:
        "Subsetting annotation GTF on {wildcards.subset}"
    threads:
        7
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 8192),
        time_min = lambda wildcards, attempt: attempt * 45
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/subset_gtf_reference/{subset}/less0.log",
        grep1 = "logs/subset_gtf_reference/{subset}/grep1.log",
        awk2 = "logs/subset_gtf_reference/{subset}/awk2.log",
        sed3 = "logs/subset_gtf_reference/{subset}/sed3.log",
        extract4 = "logs/subset_gtf_reference/{subset}/extract4.log",
        cut5 = "logs/subset_gtf_reference/{subset}/cut5.log",
        sed6 = "logs/subset_gtf_reference/{subset}/sed6.log"
    shell:
        """less {input} 2> {log.less0} | grep -v ^# 2> {log.grep1} | awk '{{print $10,$0}}' 2> {log.awk2} | sed 's/"//;s/"//;s/;/\t/' 2> {log.sed3} | {input.script} {input.subset} 1 - 1 2> {log.extract4} | cut -f2- 2> {log.cut5} |sed 's/^\s//' > {output} 2> {log.sed6}"""


"""
This rule covers lines 41 and 42 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule gffread_extraction:
    input:
        unpack(get_subread_pairs)
    output:
        temp("gffread/{database}.gffread.fa")
    message:
        "Extrating sequences on {wildcards.database} with GFFRead"
    group:
        "gffread"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 2048, 8192),
        time_min = lambda wildcards, attempt: attempt * 45
    conda:
        "../../envs/gffread.yaml"
    log:
        "gffreads/{database}.log"
    shell:
        "gffread {input.gtf} -g {input.fasta} -w {output} > {log} 2>&1"


"""
This rule covers line 43 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule concatenate_gffread_extractions:
    input:
        expand(
            "gffread/{database}.gffread.fa",
            database = ["NONCODEv5_human_hg38_lncRNA", "gencode.v27"]
        )
    output:
        temp("gffread/GENCODE_NONCODE.transcripts.tmp")
    message:
        "Concatenation of GFFread results"
    group:
        "gffread"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 256, 512),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        "gffread/concatenate.log"
    shell:
        "cat {input} > {output} 2> {log}"


"""
This rule downloads additional files that are not scripts at all
"""
rule download_resources:
    input:
        lambda wildcards: HTTP.remote(resources[wildcards.resource])
    output:
        "resources/git/{resource}"
    message:
        "Downloading {wildcards.resource} from github"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 2048),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        "logs/download/{resource}.log"
    shell:
        "wget {input} 2> {log}"


"""
This rule extracts txids from a precise salmon file.
It covers line 45 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule extract_txids:
    input:
        "../output/Salmon/TCGA-AX-A3FS-01A-11R-A22K-07_simul/quant.sf"
    output:
        temp("resources/gencode.v27.annotation.gtf.txids")
    message:
        "Extracting txids"
    threads:
        5
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 2048),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/extract_txids/less0.log",
        cut1 = "logs/extract_txids/cut1.log",
        sed2 = "logs/extract_txids/sed2.log",
        awk3 = "logs/extract_txids/awk3.log",
        sed4 = "logs/extract_txids/sed4.log"
    shell:
        """less {input} 2> {log.less0} | cut -f1 2> {log.cut1} | sed 1d 2> {log.sed2} | awk '{{print $1"\t"$1}}' 2> {log.awk3} | sed 's/|ENSG[[:graph:]]\+//' > {output} 2> {log.sed4}"""


"""
This rule covers line 46 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule gencode_subset_gene:
    input:
        "resources/gencode.v27.annotation.gtf"
    output:
        temp("resources/gencode.v27.annotation.gtf.tmp")
    message:
        "Subsetting gencode for further processes"
    threads:
        5
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 1024, 2048),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/gencode_subset_gene/less0.log",
        grep1 = "logs/gencode_subset_gene/grep1.log",
        awk2 = "logs/gencode_subset_gene/awk2.log",
        awk3 = "logs/gencode_subset_gene/awk3.log",
        sed4 = "logs/gencode_subset_gene/sed4.log"
    shell:
        """less {input} 2> {log.less0} | grep -v ^# 2> {log.grep1} | awk '$3!="gene"' 2> {log.awk2} | awk '{{print $12,$0}}' 2> {log.awk3} | sed 's/^"//;s/";\s\+/\t/' > {output} 2> {log.sed4}"""


"""
This rule covers line 47 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule match_id:
    input:
        script = "scripts/match.pl",
        gtf_tmp = "resources/gencode.v27.annotation.gtf.tmp",
        tids = "resources/gencode.v27.annotation.gtf.txids"
    output:
        temp("resources/gencode.v27.annotation.gtf.tmp1")
    message:
        "Matching extracted transcript ids with gencode GTF"
    threads:
        4
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 45
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/match_id/less0.log",
        match1 = "logs/match_id/match1.log",
        paste2 = "logs/match_id/paste2.log",
        cut3 = "logs/match_id/cut3.log"
    shell:
        "less {input.gtf_tmp} 2> {log.less0} | "
        "{input.script} {input.tids} 1 - 1 2> {log.match1} | "
        "paste - {input.gtf_tmp} 2> {log.paste2} | "
        "cut -f2,4- > {output} 2> {log.cut3}"


"""
This rule covers line 48 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule noncode_subset:
    input:
        "references/NONCODEv5_human_hg38_lncRNA.gtf"
    output:
        temp("resources/NONCODEv5_human_hg38_lncRNA.gtf.tmp")
    message:
        "Subsseting NONCODEv5_human_hg38_lncRNA"
    threads:
        4
    resources:
        mem_mb = lambda wildcards, attempt: min(512 * attempt, 1024),
        time_min = lambda wildcards, attempt: attempt * 30
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/noncode_subset/less0.log",
        grep1 = "logs/noncode_subset/grep1.log",
        awk2 = "logs/noncode_subset/awk2.log",
        sed3 = "logs/noncode_subset/sed3.log"
    shell:
        """less {input} 2> {log.less0} | grep -v ^# 2> {log.grep1} | awk '{{print $12,$0}}' 2> {log.awk2} | sed 's/^"//;s/";\s\+/\t/' > {output} 2> {log.sed3}"""


"""
This rule covers line 49 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule merge_annotations:
    input:
        noncode = "resources/NONCODEv5_human_hg38_lncRNA.gtf.tmp",
        gencode = "resources/gencode.v27.annotation.gtf.tmp1",
        script = "scripts/Exclude.pl",
        clusters = "resources/git/duplicate_clusters.tsv"
    output:
        "resources/GENCODE_NONCODE.annotation.gtf"
    message:
        "Merging annotations GencodeV27 and NONCODEv5"
    threads:
        4
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    log:
        cat0 = "logs/merge_annotations/cat0.log",
        exclude1 = "logs/merge_annotations/exclude1.log",
        cut2 = "logs/merge_annotations/cut2.log",
        awk3 = "logs/merge_annotations/awk3.log"
    shell:
        """cat {input.gencode} {input.noncode} 2> {log.cat0} | {input.script} {input.cluster} 2 - 1 2> {log.exclude1} | cut -f2- 2> {log.cut2} | awk '$7!="." && $1!~/_/' > {output} 2> {log.awk3}"""


"""
This rule covers line 51 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule get_transcripts_id_from_merged_annotations:
    input:
        gencode = "resources/gencode.v27.annotation.gtf.tmp1",
        noncode = "resources/NONCODEv5_human_hg38_lncRNA.gtf.tmp",
        script = "scripts/Exclude.pl",
        clusters = "resources/duplicate_clusters.tsv"
    output:
        temp("resources/GENCODE_NONCODE.transcripts.ids")
    message:
        "Gathering transcript identifiers from both GencodeV27 and NONCODEv5"
    threads:
        4
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    log:
        cat0 = "logs/get_transcripts_merged_annotations/cat0.log",
        exclude1 = "logs/get_transcripts_merged_annotations/exclude1.log",
        cut3 = "logs/get_transcripts_merged_annotations/cut3.log",
        awk2 = "logs/get_transcripts_merged_annotations/awk2.log"
    shell:
        """cat {input.gencode} {input.noncode} 2> {log.cat0} | {input.script} {input.cluster} 2 - 1 2> {log.exclude1} |  awk '$8!="." && $2!~/_/ && $4=="transcript"' 2> {log.awk2} | cut -f1 > {output} 2> {log.cut3}"""


"""
This rule covers line 52 of:
https://github.com/gevaertlab/RNASeq_pipeline/blob/master/prepare_reference/run_reference.sh
"""
rule merge_fasta:
    input:
        gffread = "gffread/GENCODE_NONCODE.transcripts.tmp",
        txids = "resources/GENCODE_NONCODE.transcripts.ids",
        script = "scripts/Extract.pl"
    output:
        "resources/GENCODE_NONCODE.transcripts.fa"
    message:
        "Gathering transcript identifiers from both GencodeV27 and NONCODEv5"
    threads:
        5
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    log:
        less0 = "logs/merge_fasta/less0.log",
        awk1 = "logs/merge_fasta/awk1.log",
        sed2 = "logs/merge_fasta/sed2.log",
        extract3 = "logs/merge_fasta/extract3.log",
        sed4 = "logs/merge_fasta/sed4.log"
    shell:
        """less {input.gffread} 2> {log.less0} |  awk '!/^>/ {{ printf "%s", $0; n = "\n" }} /^>/ {{ print n $0; n = "" }} END {{ printf "%s", n }}' 2> {log.awk1} | sed 'N;s/\n/\t/g;s/>/>\t/' 2> {log.sed2} | {input.script}  {input.txids} 1 - 2 2> {log.extract3} | sed 's/>\t/>/;s/\t/\n/' > {output} 2> {log.sed4}"""




################
### CLEANING ###
################
"""
This rule copies from cold to hot storage
"""
rule copy_fastq:
    input:
        lambda wildcards: copy_dict[wildcards.file]
    output:
        temp("raw_data/{file}")
    message:
        "Copying {wildcards.file} for further process"
    group:
        "Raw_Data_Cleaning"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 256, 512)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 180, 360)
        )
    wildcard_constraints:
        file = r'[^/]+'
    params:
        extra = "--verbose",
        cold_storage = ["/mnt/isilon"]
    log:
        "logs/copy/{file}.log"
    wrapper:
        f"{git}/bio/cp"


"""
This rule cleans fastq files before any other process
"""
rule trim_galore_paired:
    input:
        "raw_data/{sample}_1.fastq.gz",
        "raw_data/{sample}_2.fastq.gz"
    output:
        temp("trim_galore/fastq/{sample}_1_val_1.fq.gz"),
        "trim_galore/fastq/{sample}_1.fastq.gz_trimming_report.txt",
        temp("trim_galore/fastq/{sample}_2_val_2.fq.gz"),
        "trim_galore/fastq/{sample}_2.fastq.gz_trimming_report.txt"
    message:
        "Cleaning (pair ended) {wildcards.sample} with Trim_Galore"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(4096 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 480
    group:
        "Raw_Data_Cleaning"
    conda:
        "../../envs/trim_galore.yaml"
    params:
        extra = config["params"].get(
            "trim_galore_extra",
            "-q 20 --stringency 3 --gzip --length 20"
        )
    log:
        "logs/trim_galore/{sample}.log"
    shell:
        " trim_galore "
        " {params.extra} "
        " --paired "
        " -o trim_galore/fastq/ "
        " {input} "
        " > {log} 2>&1 "


rule trim_galore_single:
    input:
        "raw_data/{sample}_1.fastq.gz"
    output:
        temp("trim_galore/fastq/{sample}_1_trimmed.fq.gz"),
        "trim_galore/fastq/{sample}_1.fastq.gz_trimming_report.txt"
    message:
        "Cleaning (single ended) {wildcards.sample} with Trim_Galore"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(4096 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 480
    group:
        "Raw_Data_Cleaning"
    conda:
        "../../envs/trim_galore.yaml"
    params:
        extra = config["params"].get(
            "trim_galore_extra",
            "-q 20 --stringency 3 --gzip --length 20"
        )
    log:
        "logs/trim_galore/{sample}.log"
    shell:
        " trim_galore "
        " {params.extra} "
        " -o trim_galore/fastq/ "
        " {input} "
        " > {log} 2>&1 "


################
### KALLISTO ###
################

"""
This rule indexes the fasta sequence used by Kallisto

Kallisto version is provided here: https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh#L2
Default arguments were used: https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/prepare_reference/run_reference.sh#L58
"""
rule kallisto_index:
    input:
        fasta = fasta_path
    output:
        index = "kallisto/index/gencode.v27.transcripts.idx"
    message:
        "Building index from {input.fasta} with Kallisto"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(5120 * attempt, 10240),
        time_min = lambda wildcards, attempt: attempt * 240
    conda:
        "../../envs/kallisto.yaml"
    params:
        extra = config['params'].get(
            "kallisto_index_extra",
            ""
        )
    priority:
        50
    log:
        "logs/kallisto/index.log"
    shell:
        """
	mkdir -p kallisto/index
	kallisto index --index={output.index} {params.extra} {input.fasta} >{log} 2>&1
	"""

"""
This rule estimates fragment size for single ended samples
"""
rule fragment_size_estimate:
    input:
        "trim_galore/fastq/{sample}_1_trimmed.fq.gz"
    output:
        temp("fragment_size/{sample}.txt")
    message:
        "Estimating fragment size since Kallisto does not want to do it for "
        "single end. Performed on {wildcards.sample}"
    threads:
        2
    resources:
        time_min = 15,
        mem_mb = 512
    params:
        nb_reads = config["params"].get("nb_read_fragment_size_estimate", 100)
    log:
        "logs/fragment_size/{sample}.log"
    shell:
        """awk 'BEGIN {{FS="\\t"; READ_LEN_SUM=0; READ_SUM_SQUARE=0}} NR % 4 == 0 {{READ_LEN_SUM=READ_LEN_SUM+length($0); READ_SUM_SQUARE=READ_SUM_SQUARE+(length($0)^2)}} END {{READS_NB=NR/4; MEAN=READ_LEN_SUM/READS_NB; SD=sqrt(READ_SUM_SQUARE/READS_NB-MEAN^2); print "-l "MEAN" -s "SD}}' <(gunzip -c {input}) > {output} 2> {log}"""


"""
This rule quantifies fastq reads with Kallisto

Kallisto version is provided here: https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/run_workflow.sh#L2
Kallisto quant arguments for unstranded library is provided here: https://github.com/gevaertlab/RNASeq_pipeline/blob/32b740603a0bba11dafb4b1facc1b2a92dd030dc/run_workflow.sh#L78
"""
rule kallisto_quant:
    input:
        unpack(get_input_fq_files_kallisto),
        index = "kallisto/index/gencode.v27.transcripts.idx",
        gtf = gtf_path
    output:
        expand(
            "kallisto/quant/{sample}/{content}",
            content = ["abundance.h5", "abundance.tsv", "run_info.json"],
            allow_missing = True
        )
    message:
        "Quantifying {wildcards.sample} with Kallisto"
    threads:
        min(config.get("threads", 1), 20)
    priority:
        50
    resources:
        mem_mb = lambda wildcards, attempt: min(8192 * attempt, 16384),
        time_min = lambda wildcards, attempt: attempt * 480
    conda:
        "../../envs/kallisto.yaml"
    params:
        extra = config["params"].get(
            "kallisto_quant_extra",
            ""
        ),
        single = lambda w: "" if is_paired(w) else "--single",
        fragment_size_estimate = (
            lambda w: (
                "" if is_paired(w) else f"$(cat fragment_size/{w.sample}.txt)"
            )
        )
    log:
        "logs/kallisto/quant/{sample}.log"
    shell:
        "mkdir -p kallisto/quant/{wildcards.sample} && "
        "kallisto quant {params.extra} "
        "--gtf={input.gtf} --threads={threads} "
        "--output-dir=kallisto/quant/{wildcards.sample} "
        "{params.single} {params.fragment_size_estimate} "
        "--index={input.index} {input.fastq} >{log} 2>&1"


"""
This rule archives the quantifications in order to gain space
"""
rule archive_quantification:
    input:
        quant = expand(
            "kallisto/quant/{sample}/{content}",
            content = ["abundance.h5", "abundance.tsv", "run_info.json"],
            allow_missing = True
        )
    output:
        archive = report(
            "kallisto/quant/{sample}.tar.bz2",
            caption="../../reports/kallisto_archive.rst",
            category="Kallisto"
        )
    message:
        "Compressing kallisto quantification result for {wildcards.sample}"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(2048 * attempt, 8192),
        time_min = lambda wildcards, attempt: attempt * 60
    conda:
        "../../envs/bash_perl.yaml"
    params:
        extra = "-cvjf"
    log:
        "logs/kallisto/archive/{sample}.log"
    shell:
        " tar "
        " {params.extra} "
        " {output.archive} "
        " {input.quant} "
        " 2> {log}"

################
### TXIMPORT ###
################
"""
This rule computes tx2gene, since original pipeline does not provide Any
"""
rule tx2gene:
    input:
        gtf = gtf_path
    output:
        tsv = temp("reference/tx2gene.tsv")
    message:
        "Building table of correspondacy between genes and transcript names"
    group:
        "tx2gene"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(2048 * attempt, 10240),
        time_min = lambda wildcards, attempt: attempt * 60
    log:
        "logs/tx2gene/tr2gene.log"
    wrapper:
        f"{git}/bio/gtf/tx2gene"


"""
This rule cuts the values to what tximport requires
"""
rule tximport_tx2gene:
    input:
        "reference/tx2gene.tsv"
    output:
        "tximport/tx2gene.tsv"
    message:
        "Cutting correspondacy table to fit tximport requirements"
    group:
        "tx2gene"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(1024 * attempt, 10240),
        time_min = lambda wildcards, attempt: attempt * 30
    log:
        "logs/tximport/tx2gene.log"
    conda:
        "../../envs/bash_perl.yaml"
    shell:
        """awk '{{print $2"\t"$1}}' {input} > {output} 2> {log}"""


"""
This rule imports counts in R with tximport
"""
rule tximport:
    input:
        quant = expand(
            "kallisto/quant/{sample}/abundance.h5",
            sample=design["Sample_id"].tolist()
        ),
        tx2gene = "tximport/tx2gene.tsv"
    output:
        txi = report(
            "tximport/txi.rds",
            caption="../../reports/tximport.rst",
            category="tximport"
        ),
        txt = report(
            "tximport/quantification.csv",
            caption="../../reports/tximport_quant.rst",
            category="tximport"
        )
    message:
        "Importing counts in R with tximport"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(attempt * 20480, 81920),
        time_min = lambda wildcards, attempt: attempt * 180
    conda:
        "../../envs/tximport.yaml"
    log:
        "logs/tximport/txi.log"
    script:
        "../../scripts/tximport.R"
