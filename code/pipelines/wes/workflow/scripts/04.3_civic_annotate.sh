#!/bin/bash

while getopts ":a:b:c:d:e:f:m:n:o:t:l:" opt; do
    case $opt in
        a) table_alt+=("$OPTARG")
	    ;;
        b) table_cln="$OPTARG"
            ;;
        c) table_gen="$OPTARG"
            ;;
        d) table_pre="$OPTARG"
            ;;
        e) table_run="$OPTARG"
            ;;
        f) table_pos="$OPTARG"
            ;;
        m) code_dir="$OPTARG"
            ;;
        n) civic="$OPTARG"
            ;;
        o) rules="$OPTARG"
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

# preprocess table of alterations before annotating with CIViC
printf -- "-INFO: preparing table before annotating with CIViC...\n"
python -u workflow/scripts/04.1_civic_preprocess.py \
    --table_alt ${table_alt} \
    --table_cln ${table_cln} \
    --table_gen ${table_gen} \
    --gen_gene_name name \
    --category ${category} \
    --output ${table_pre}

# run CIViC annotator
printf -- "-INFO: running CIViC annotator...\n"
python -u ${code_dir}/civic.py \
    --input ${table_pre} \
    --civic ${civic} \
    --rules ${rules} \
    --category ${category} \
    --output ${table_run}

# postprocess annotations
printf -- "-INFO: postprocess CIViC annotations...\n"
python -u workflow/scripts/04.2_civic_postprocess.py \
    --input ${table_run} \
    --output ${table_pos}
