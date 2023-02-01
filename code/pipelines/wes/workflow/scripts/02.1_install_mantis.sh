#!/bin/bash

while getopts ":c:f:s:n:b:v:d:" opt; do
    case $opt in
        c) code_dir="$OPTARG"
            ;;
        f) fasta="$OPTARG"
            ;;
        s) ms_dir="$OPTARG"
            ;;
        n) targets_name+=("$OPTARG")
            ;;
        b) targets_bed+=("$OPTARG")
            ;;
        v) buildver="$OPTARG"
            ;;
        d) ms_bed="$OPTARG"
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
printf "Hello! This script is intended to build lists of microsatellites using a reference genomes and one or more target files\n"

#### Mantis code ####

local_file="${code_dir}/mantis.py"

if [[ -f ${local_file} ]]
then
    printf "* the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    url="https://github.com/OSU-SRLab/MANTIS"
    printf "* downloading repository from %s...\n" "$url"

    git clone ${url} ${code_dir}
fi

# Build microsatellites list from reference genome
ms_genome_unfixed="${ms_dir}/${buildver}_microsatellites_unfixed.bed"
ms_genome="${ms_dir}/${buildver}_microsatellites.bed"

if [[ -f ${ms_genome} ]]
then
    printf "* the output file %s already exists! Delete it if you want to overwrite it\n" "$ms_genome"
else
    printf "* making RepeatFinder ...\n"

    cwd=$(pwd)
    cd ${code_dir}/tools
    make RepeatFinder
    cd ${cwd}

    printf "* running RepeatFinder on %s...\n" "${fasta}"
    mkdir -p ${ms_dir}

    echo "* ./${code_dir}/tools/RepeatFinder –i ${fasta} -o ${ms_genome_unfixed}"
    ./${code_dir}/tools/RepeatFinder -i ${fasta} -o ${ms_genome_unfixed}

    echo "* python ./${code_dir}/tools/fix_RF_bed_output.py -i ${ms_genome_unfixed} -o ${ms_genome}"
    python ./${code_dir}/tools/fix_RF_bed_output.py -i ${ms_genome_unfixed} -o ${ms_genome}

    rm ${ms_genome_unfixed}
fi

# Intersect microsatellites reference genome list with some list of microsatellites
ms_name=$(basename ${ms_bed})
ms_genome_bed="${ms_dir}/${buildver}_microsatellites_${ms_name}"

printf "* intersecting %s with  %s...\n" "${ms_genome}" "${ms_bed}"
echo "* bedtools intersect -a ${ms_genome} -b ${ms_bed} -wa > ${ms_genome_bed}"

bedtools intersect -a ${ms_genome} -b ${ms_bed} -wa > ${ms_genome_bed}

# fix start-1 error
awk -F $"\t" 'BEGIN {OFS=FS} {$2=$2+1; print}' ${ms_genome_bed} > a
mv a ${ms_genome_bed}

# Intersect microsatellites reference genome list with target lists

for ((i=0;i<${#targets_name[@]};i++))
do
    target_bed=${targets_bed[$i]}
    target_name=${targets_name[$i]}
    target_ms_genome="${ms_dir}/${buildver}_microsatellites_${target_name}.bed"

    printf "* intersecting %s with  %s...\n" "${ms_genome}" "${target_bed}"
    echo "* bedtools intersect -a ${ms_genome} -b ${target_bed} -wa > ${target_ms_genome}"

    bedtools intersect -a ${ms_genome} -b ${target_bed} -wa > ${target_ms_genome}

    # fix start-1 error
    awk -F $"\t" 'BEGIN {OFS=FS} {$2=$2+1; print}' ${target_ms_genome} > a
    mv a ${target_ms_genome}
done

printf "All done!\n"
