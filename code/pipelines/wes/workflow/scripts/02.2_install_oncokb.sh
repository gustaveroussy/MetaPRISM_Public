#!/bin/bash

while getopts ":c:d:g:" opt; do
    case $opt in
        c) code_dir="$OPTARG"
            ;;
        d) data_dir="$OPTARG"
            ;;
        g) gene_list="$OPTARG"
            ;;
        \?) echo "Invalid option -$OPTARG" >&2
            exit 1
            ;;
    esac

    case $OPTARG in
        -*) echo "Option $opt needs a valid argument"
            exit 1
            ;;
    esac
done

module load gcc

# message
printf "Hello! This script is intended to download oncokb annotator code from github and some tables from the oncokb website.\n"

#### Onockb code ####
local_file="${code_dir}/MafAnnotator.py"

if [[ -f ${local_file} ]]
then
    printf "* the output file %s already exists! Delete it if you want to overwrite it\n" "${local_file}"
else
    if [[ -d ${code_dir} ]]
    then
        rm -r ${code_dir}
    fi

    url="https://github.com/oncokb/oncokb-annotator"
    printf "* downloading repository from %s...\n" "$url"

    git clone ${url} ${code_dir}
fi

#### Download tables ####

local_file_1="${data_dir}/cancerGeneList.tsv"

if [[ -f ${local_file_1} ]]
then
    printf "* the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file_1"
else
    if [[ ! -f ${gene_list} ]]
    then
        url="https://www.oncokb.org/cancerGenes"
        printf "* the output file %s does not exist! Download it manually from %s\n" "$local_file_1" "$url"
        exit 1
    else
        cp ${gene_list} ${local_file_1}
    fi
fi

local_file_2="${data_dir}/cancerGeneList_oncokb_annotated.tsv"

if [[ -f ${local_file_2} ]]
then
    printf "* the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file_2"
else
    printf "* selecting Oncokb annotated genes %s\n"
    awk -F '\t' 'NR==1; NR>1{ if ($8=="Yes") {print}}' ${local_file_1} > ${local_file_2}
fi
