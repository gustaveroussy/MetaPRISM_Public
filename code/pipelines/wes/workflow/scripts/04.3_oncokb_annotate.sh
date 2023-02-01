#!/bin/bash

while getopts ":a:b:c:d:g:e:f:m:o:k:t:l:" opt; do
    case $opt in
        a) table_alt+=("$OPTARG")
	    ;;
        b) table_cln="$OPTARG"
            ;;
        c) table_gen="$OPTARG"
            ;;
        d) table_alt_pre="$OPTARG"
            ;;
        g) table_cln_pre="$OPTARG"
            ;;
        e) table_run="$OPTARG"
            ;;
        f) table_pos="$OPTARG"
            ;;
        m) code_dir="$OPTARG"
            ;;
        o) rules="$OPTARG"
            ;;
        k) token="$OPTARG"
            ;;
        t) category="$OPTARG"
            ;;
        l) log="$OPTARG"
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

exec 3>&1 4>&2 >>${log} 2>&1

# convert array to space-separated values
table_alt=$(printf "%s" "${table_alt[*]}")

# preprocess table of alterations before annotating with OncoKB
printf -- "-INFO: preparing table before annotating with OncoKB...\n"
python -u workflow/scripts/04.1_oncokb_preprocess.py \
    --table_alt ${table_alt} \
    --table_cln ${table_cln} \
    --table_gen ${table_gen} \
    --gen_gene_name "Hugo Symbol" \
    --category ${category} \
    --output_alt ${table_alt_pre} \
    --output_cln ${table_cln_pre}

# run OncoKB-annotator
printf -- "-INFO: running OncoKB-annotator...\n"

if [[ ${category} == "cna" ]]; then
    python -u ${code_dir}/CnaAnnotator.py \
        -i ${table_alt_pre} \
        -c ${table_cln_pre} \
        -b ${token} \
        -o ${table_run}
elif [[ ${category} == "mut" ]]; then
    python -u ${code_dir}/MafAnnotator.py \
        -i ${table_alt_pre} \
        -c ${table_cln_pre} \
        -b ${token} \
        -q Genomic_Change \
        -o ${table_run}
else
    printf -- "-ERROR: unrecognized value of --category '%s'\n" "${category}"
    exit -1
fi

# postprocess annotations
printf -- "-INFO: postprocess OncoKB annotations...\n"
python -u workflow/scripts/04.2_oncokb_postprocess.py \
    --input ${table_run} \
    --rules ${rules} \
    --category ${category} \
    --output ${table_pos}
