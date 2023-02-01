#!/bin/bash

while getopts ":b:x:y:" opt; do
    case $opt in
	b) bam="$OPTARG"
	    ;;
	x) fastq_1="$OPTARG"
	    ;;
	y) fastq_2="$OPTARG"
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


bam_filtered=${bam%.bam}.filtered.bam
bam_sorted=${bam%.bam}.filtered.sorted.bam

echo "filtering reads..."
echo "samtools view -h -b -f 3 -F 3840 ${bam} > ${bam_filtered}"
samtools view -h -b -f 3 -F 3840 ${bam} > ${bam_filtered}

echo "sorting reads..."
echo "samtools sort -n -o ${bam_sorted} ${bam_filtered}"
samtools sort -n -o ${bam_sorted} ${bam_filtered}

parent_1=$(dirname $fastq_1)
if [[ ! -d "${parent_1}" ]]; then
    mkdir -p ${parent_1}
fi

parent_2=$(dirname $fastq_2)
if [[ ! -d "${parent_2}" ]]; then
    mkdir -p ${parent_2}
fi

echo "converting bam to fastq..."
echo "bedtools bamtofastq -i ${bam_sorted} -fq ${fastq_1} -fq2 ${fastq_2}"
bedtools bamtofastq -i ${bam_sorted} -fq ${fastq_1} -fq2 ${fastq_2}

echo "remove ${bam_filtered} and ${bam_sorted}..."
rm ${bam_filtered}
rm ${bam_sorted}
