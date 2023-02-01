#!/bin/bash

while getopts ":d:i:n:c:t:a:" opt; do
    case $opt in
	d) dir="$OPTARG"
	    ;;
	i) bam_id="$OPTARG"
	    ;;
	n) bam_name="$OPTARG"
	    ;;
	c) gdc_client="$OPTARG"
	    ;;
	t) gdc_token="$OPTARG"
	    ;;
	a) threads="$OPTARG"
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

raw_bam=${bam_id}/${bam_name}
raw_bam_bai=${bam_id}/${bam_name}.bai
raw_bam_logs=${bam_id}/logs

echo "downloading BAM file to ${raw_bam} using gdc-client..."
./${gdc_client} download -t ${gdc_token} -n ${threads} ${bam_id}

rm ${raw_bam_bai}
echo "removed file ${raw_bam_bai}"

rm -r ${raw_bam_logs}
echo "removed folder ${raw_bam_logs}"

if [[ ! -d "${dir}" ]];then
    mkdir -p ${dir}
    echo "created directory ${dir}"
fi

new_bam=${dir}/${bam_name}
mv ${raw_bam} ${new_bam}
echo "moved bam file to ${new_bam}"

rm -r ${bam_id}
echo "removed folder ${bam_id}"
