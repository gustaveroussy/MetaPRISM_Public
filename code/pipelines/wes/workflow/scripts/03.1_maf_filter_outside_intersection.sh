#!/bin/bash

while getopts ":d:v:t:n:b:o:" opt; do
    case $opt in
        d) dir="$OPTARG"
            ;;
        v) vcf="$OPTARG"
            ;;
        t) tsample="$OPTARG"
            ;;
        n) nsample="$OPTARG"
            ;;
        b) bed="$OPTARG"
            ;;
        o) output="$OPTARG"
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

header="${dir}/${tsample}_vs_${nsample}_header.vcf"
inbed="${dir}/${tsample}_vs_${nsample}_content_inbed.vcf"
offbed="${dir}/${tsample}_vs_${nsample}_content_offbed.vcf"
offbed_ann="${dir}/${tsample}_vs_${nsample}_content_offbed_ann.vcf"

bedtools intersect -a ${vcf} -b ${bed} -wa | sort | uniq > ${inbed}
bedtools intersect -a ${vcf} -b ${bed} -v | sort | uniq > ${offbed}
zgrep "^#" ${vcf} > ${header}
gawk 'BEGIN{FS="\t"; OFS="\t"} {if ($7=="PASS") $7="OFF_TARGETS_INTERSECTION"; else $7="OFF_TARGETS_INTERSECTION;"$7; print}' ${offbed} > ${offbed_ann}
sed -i '3 i ##FILTER=<ID=OFF_TARGETS_INTERSECTION,Description='"\"Variant is outside the bed file ${bed}.\""'>' ${header}
cat ${header} ${inbed} ${offbed_ann} | bgzip -c > ${output}
rm ${header} ${inbed} ${offbed} ${offbed_ann}
picard SortVcf I=${output} O=${output}
