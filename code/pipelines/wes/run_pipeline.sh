#!/bin/bash

usage() { echo "$0 Usage:" && grep " .)\ #" $0; exit 0; }
one_run="yes"
aggregate="no"
dry_run="no"
jobs=100
r_folder=results

while getopts ":a:b:i:j:r::ogn h" opt; do
  case $opt in
    a) # Minimum batch index to be run. Use -o option if you want to run the pipeline in one batch. If running in multiple batches, you must first split the file config/samples.tsv into as many files `config/samples_${index}` as the number of batches. Ignored if -i is used or if -o is used. 
      batch_min="$OPTARG"
      ;;
    b) # Maximum batch index to be run. Ignored if -i is used or if -o is used.
      batch_max="$OPTARG"
      ;;
    i) # Batch index to be run. Ignored if -o is used.
      batch_idx="$OPTARG"
      ;;
    j) # Default=100. Maximum number of jobs that may be run in paralle by snakemake.
      jobs="$OPTARG"
      ;;
    r) # Default=results. Relative or absolute path to results folder. This should mirror the value in rules/common.smk.
      r_folder="$OPTARG"
      ;;
    o) # If -o option is used, the pipeline is run in one batch using all samples in config/samples.tsv.
      one_run="yes"
      ;;
    g) # If -g option is used, the aggregation rules are run locally after the last bacth has been processed. Ignored if -o is used.
      aggregate="yes"
      ;;
    n) # If -n option is used, dry-run.
      dry_run="yes"
      ;;
    h) # Display help.
      usage
      ;;
    *) echo "Invalid option -$OPTARG" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ ${one_run} == "yes" ]]; then
  if [[ ${dry_run} == "yes" ]]; then
      snakemake -s workflow/Snakefile --profile slurm --jobs ${jobs} --rerun-incomplete -n
  else
      snakemake -s workflow/Snakefile --profile slurm --jobs ${jobs} --rerun-incomplete
  fi
else
  if [[ -z "${batch_idx}" ]]; then
    batch_indices=("${batch_idx}")
  else
    batch_indices=($(seq ${batch_min} 1 ${batch_max}))
  fi

  config_samples=config/samples.tsv
  config_samples_bu=config/samples.bu.tsv
  config_yaml=config/config.yaml
  config_yaml_bu=config/config.bu.yaml
  
  if [[ -f "${config_samples}" ]]; then
    printf "INFO: backing up the file %s as %s\n" "${config_samples}" "${config_samples_bu}"
    mv ${config_samples} ${config_samples_bu}
  fi

  if [[ -f "${config_yaml}" ]]; then
    mv ${config_yaml} ${config_yaml_bu}
  fi

  if ! grep -R '^skip_aggregate:' ${config_yaml} > /dev/null ; then
    printf "skip_aggregate: true" >> ${config_yaml} 
  else
    sed -i -e "s/skip_aggregate:.*/skip_aggregate: true/" ${config_yaml} 
  fi

  for batch_index in "${batch_indices[@]}"
  do
    printf "\n\nRunning batch ${batch_index}/[${batch_min}, ${batch_max}] ...\n\n"
    cp config/samples_${batch_index}.tsv ${config_yaml}
    if [[ ${dry_run} == "yes" ]]; then
      snakemake -s workflow/Snakefile --profile slurm --jobs ${jobs} --rerun-incomplete -n
    else
      snakemake -s workflow/Snakefile --profile slurm --jobs ${jobs} --rerun-incomplete
    fi
  done

  printf "\nAll batches have been processed! Restoring backup files.\n"

  if [[ -f "${config_samples}" ]]; then
    rm ${config_samples}
  fi

  if [[ -f "${config_samples_bu}" ]]; then
    mv ${config_samples_bu} ${config_samples}
  fi

  if [[ -f "${config_yaml_bu}" ]]; then
    mv ${config_yaml_bu} ${config_yaml}
  fi

  if [[ ${aggregate} == "yes" ]]; then
    printf "Running aggregation rules locally...\n"
    printf "WARNING 1: this may be inappropriate if you are currently running this script from the login node of a shared HPC and the aggregations are computationally-intensive.\n"
    printf "WARNING 2: the local aggregation will aggregate all the files found in each results folder to be aggregated. It may not correspond to your current selection of samples.\n"

    # somatic_maf_filters_aggregate
    folder=${r_folder}/calling/somatic_maf_filters
    if [[ -d "${folder}" ]]; then
      output=${r_folder}/aggregate/somatic_maf/somatic_calls_filters_prefinal.tsv.gz
      mkdir -p  ${r_folder}/aggregate/somatic_maf
      if [[ ${dry_run} == "no" ]]; then
        printf "aggregating contents of %s to %s...\n" "${folder}" "${output}"
        zcat ${folder}/*.tsv.gz | sed -n '1p;/^CHROM/ !p' | gzip > ${output}
      else
        printf "(dry-run) aggregating contents of %s to %s...\n" "${folder}" "${output}"
      fi
    fi


    # somatic_maf_aggregate
    folder=${r_folder}/annotation/somatic_maf
    if [[ -d "${folder}" ]]; then
      output=${r_folder}/aggregate/somatic_maf/somatic_calls_prefinal.tsv.gz
      mkdir -p  ${r_folder}/aggregate/somatic_maf

      if [[ ${dry_run} == "no" ]]; then
        printf "aggregating contents of %s to %s...\n" "${folder}" "${output}"
        files=( $folder/*.maf )
        if [[ ${#files[@]} == 1 ]]; then
            cat ${files[0]} | gzip > ${output}
        else
            { cat ${files[@]:0:1}; grep -v "^##\|Tumor_Sample_Barcode" ${files[@]:1}; } | gzip > ${output}
        fi
      else
        printf "(dry-run) aggregating contents of %s to %s...\n" "${folder}" "${output}"
      fi
    fi

    # somatic_maf_oncokb_aggregate
    folder=${r_folder}/annotation/somatic_maf_oncokb
    if [[ -d "${folder}" ]]; then
      output=${r_folder}/aggregate/somatic_maf/somatic_calls_oncokb_prefinal.tsv.gz
      mkdir -p  ${r_folder}/aggregate/somatic_maf

      if [[ ${dry_run} == "no" ]]; then
        printf "aggregating contents of %s to %s...\n" "${folder}" "${output}"
        files=( $folder/*.maf )
        if [[ ${#files[@]} == 1 ]]; then
            cat ${files[0]} | gzip > ${output}
        else
            { cat ${files[@]:0:1}; grep -v "^##\|Tumor_Sample_Barcode" ${files[@]:1}; } | gzip > ${output}
        fi
      else
        printf "(dry-run) aggregating contents of %s to %s...\n" "${folder}" "${output}"
      fi
    fi

    # somatic_maf_civic_aggregate
    folder=${r_folder}/annotation/somatic_maf_civic
    if [[ -d "${folder}" ]]; then
      output=${r_folder}/aggregate/somatic_maf/somatic_calls_civic_prefinal.tsv.gz
      mkdir -p  ${r_folder}/aggregate/somatic_maf

      if [[ ${dry_run} == "no" ]]; then
        printf "aggregating contents of %s to %s...\n" "${folder}" "${output}"
        files=( $folder/*.maf )
        if [[ ${#files[@]} == 1 ]]; then
            cat ${files[0]} | gzip > ${output}
        else
            { cat ${files[@]:0:1}; grep -v "^##\|Tumor_Sample_Barcode" ${files[@]:1}; } | gzip > ${output}
        fi
      else
        printf "(dry-run) aggregating contents of %s to %s...\n" "${folder}" "${output}"
      fi
    fi
  fi
fi
