#!/bin/bash

while getopts ":c:a:u:v:o:" opt; do
    case $opt in
        c) vep_cache="$OPTARG"
            ;;
        u) cadd_version="$OPTARG"
            ;;
        v) dbnsfp_version="$OPTARG"
            ;;
        a) assembly="$OPTARG"
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

printf "Hello! This script is intended to retrieve databases required for running VEP with some plugins\n"
printf "=%.0s" {1..80}
printf  "\n\n"
printf "* databases will be be added to the vep cache folder: %s\n" "$vep_cache"
printf "* the reference genome assembly is: %s\n" "$assembly"
printf "* the dbnsfp version required is: %s\n\n" "$dbnsfp_version"

tmp_dir="/mnt/beegfs/scratch/tmp"
local_dir="${vep_cache}/Plugins"

#### CADD ####

beegfs_dir="/mnt/beegfs/database/bioinfo/Index_DB/CADD/${cadd_version}/${assembly}"
beegfs_file="/mnt/beegfs/database/bioinfo/Index_DB/CADD/${cadd_version}/${assembly}/whole_genome_SNVs.tsv.gz"
local_file="${local_dir}/whole_genome_SNVs_${cadd_version}_${assembly,,}.tsv.gz"

if [[ -f ${local_file} ]]
then
    printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    cwd=$(pwd)
    cd ${vep_cache}

    if [[ -f ${beegfs_file} ]]
    then
        printf "** using existing file %s...\n" "$beegfs_file"
        cp ${beegfs_file} ${local_file}
        cp ${beegfs_file}.tbi ${local_file}.tbi
    else
        url="https://kircherlab.bihealth.org/download/CADD/v${cadd_version}/${assembly}/whole_genome_SNVs.tsv.gz"
        printf "** downloading file from %s...\n" "$url"

        cwd=$(pwd)
        cd ${vep_cache}
        curl ${url} > ${local_file}

        printf "** building index ...\n"
        tabix -s 1 -b 2 -e 2 ${local_file}
    fi

    cd ${cwd}
    printf "** CADD file saved at %s\n" "${local_file}"
fi


beegfs_dir="/mnt/beegfs/database/bioinfo/Index_DB/CADD/${cadd_version}/${assembly}"
beegfs_file="/mnt/beegfs/database/bioinfo/Index_DB/CADD/${cadd_version}/${assembly}/InDels.tsv.gz"
local_file="${local_dir}/InDels_${cadd_version}_${assembly,,}.tsv.gz"

if [[ -f ${local_file} ]]
then
    printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    cwd=$(pwd)
    cd ${vep_cache}

    if [[ -f ${beegfs_file} ]]
    then
        printf "** using existing file %s...\n" "$beegfs_file"
        cp ${beegfs_file} ${local_file}
        cp ${beegfs_file}.tbi ${local_file}.tbi
    else
        url="https://kircherlab.bihealth.org/download/CADD/v${cadd_version}/${assembly}/InDels.tsv.gz"
        printf "** downloading file from %s...\n" "$url"

        cwd=$(pwd)
        cd ${vep_cache}
        curl ${url} > ${local_file}

        printf "** building index ...\n"
        tabix -s 1 -b 2 -e 2 ${local_file}
    fi

    cd ${cwd}
    printf "** CADD file saved at %s\n" "${local_file}"
fi


#### dbNSFP #####

beegfs_dir="/mnt/beegfs/database/bioinfo/Index_DB/dbNSFP/${dbnsfp_version}/${assembly}"
beegfs_file="${beegfs_dir}/dbNSFP${dbnsfp_version}a.txt.gz"
local_file="${local_dir}/dbNSFP${dbnsfp_version}a_${assembly,,}.gz"

printf "1. Collecting database for dbNSFP v${dbnsfp_version}a\n"
printf "+%.0s" {1..80}
printf  "\n\n"

if [[ -f ${local_file} ]]
then
    printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    if [[ -f ${beegfs_file} ]]
    then
        printf "** using existing file %s...\n" "$beegfs_file"

        zcat ${beegfs_file} | head -n1 > ${local_dir}/h
        if [ ${assembly} = "GRCh37" ]
        then
            zgrep -h -v "^#chr" ${beegfs_file} |
                awk '$8 != "." ' |
                sort -T ${tmp_dir} -k8,8 -k9,9n - |
                cat ${local_dir}/h - |
                bgzip -c > ${local_file}

            printf "** building index ...\n"
            tabix -s 8 -b 9 -e 9 ${local_file}
        else
            zgrep -h -v "^#chr" ${beegfs_file} |
                sort -T ${tmp_dir} -k1,1 -k2,2n - |
                cat ${local_dir}/h - |
                bgzip -c > ${local_file}

            printf "** building index ...\n"
            tabix -s 1 -b 2 -e 2 ${local_file}
        fi

        rm ${local_dir}/h
    else
        url="ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP${dbnsfp_version}a.zip"
        printf "** downloading file %s...\n" "$url"

        cwd=$(pwd)
        cd ${vep_cache}
        wget ${url}
        unzip dbNSFP${dbnsfp_version}a.zip
        zcat dbNSFP${dbnsfp_version}a_variant.chr1.gz | head -n1 > h

        if [ ${assembly} = "GRCh37" ]
        then
            zgrep -h -v "^#chr" dbNSFP${dbnsfp_version}a_variant.chr* |
                awk '$8 != "." ' |
                sort -T ${tmp_dir} -k8,8 -k9,9n - |
                cat h - |
                bgzip -c > dbNSFP${dbnsfp_version}a_grch37.gz 

            printf "** building index ...\n"
            tabix -s 8 -b 9 -e 9 dbNSFP${dbnsfp_version}a_grch37.gz
        else
            zgrep -h -v "^#chr" dbNSFP${dbnsfp_version}a_variant.chr* |
                sort -T ${tmp_dir} -k1,1 -k2,2n - |
                cat h -|
                bgzip -c > dbNSFP${dbnsfp_version}a_grch38.gz

            printf "** building index ...\n"
            tabix -s 1 -b 2 -e 2 dbNSFP${dbnsfp_version}a_grch38.gz
        fi

        rm h
        cd ${cwd}
    fi

    printf "** dbNSFP file saved at %s\n" "${local_file}"
    printf "** dbNSFP index file saved at %s\n" "${local_file}.tbi"
fi


#### LoFtool #####

local_file="${local_dir}/LoFtool_scores.txt"

if [[ -f ${local_file} ]]
then
    printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    url="https://github.com/Ensembl/VEP_plugins/blob/release/${release}/LoFtool_scores.txt"
    printf "** downloading file from %s...\n" "$url"

    cwd=$(pwd)
    cd ${vep_cache}
    wget ${url}
    cd ${cwd}

    printf "** LoFtool file saved at %s\n" "${local_file}"
fi

# #### PON_P2 #####
# 
# local_file="${local_dir}/ponp2.py"
# 
# if [[ -f ${local_file} ]]
# then
#     printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
# else
#     url="http://structure.bmc.lu.se/PON-P2/pluginSupportDownload/"
#     printf "** downloading file from %s...\n" "$url"
# 
#     cwd=$(pwd)
#     cd ${vep_cache}
#     curl ${url} > ${local_file}
#     cd ${cwd}
# 
#     printf "** PON_P2 file saved at %s\n" "${local_file}"
# fi

#### Condel ####

local_file="${local_dir}/Condel/config/condel_SP.conf"

if [[ -f ${local_file} ]]
then
    printf "** the output file %s already exists! Delete it if you want to overwrite it\n" "$local_file"
else
    url="https://github.com/Ensembl/VEP_plugins"
    printf "** cloning github repository %s...\n" "$url"

    cwd=$(pwd)
    cd ${vep_cache}
    git clone ${url}
    mv VEP_Plugins/config/Condel Plugins
    rm -rf VEP_Plugins
    path=$(pwd)/Plugins/config/Condel/
    sed -i 's\path/to/config/Condel/\${path}'${path}/condel_SP.conf
    cd ${cwd}

    printf "** Condel file saved at %s\n" "${local_file}"
fi

# done

printf "* all done!"
