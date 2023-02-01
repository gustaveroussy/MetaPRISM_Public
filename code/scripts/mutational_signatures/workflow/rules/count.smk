rule count_matrix_split_by_mut:
    benchmark:
        "%s/count_matrix_split_by_mut_{mmode}_{cohort}.tsv" % B_FOLDER
    input:
        mut=lambda w:config["data"]["mut"][w.cohort],
        mut_sum=lambda w:config["data"]["mut_sum"][w.cohort],
        env="../common/logs/setup_conda_signatures.done"
    conda:
        config["setup"]["Signatures"]
    log:
        "%s/count_matrix_split_by_mut_{mmode}_{cohort}.log" % L_FOLDER
    params:
        min_mut = config["count"]["mutation_min_per_tumor"]
    output:
        all_mut="%s/counts_mutations/counts_mutations_{mmode}_all_mut_{cohort}.tsv" % R_FOLDER,
        min_mut="%s/counts_mutations/counts_mutations_{mmode}_min_mut_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    threads: 1
    shell:
        """Rscript workflow/scripts/01_count_matrix.R \
            --cohort {wildcards.cohort} \
            --mut {input.mut} \
            --mut_sum {input.mut_sum} \
            --min_mut {params.min_mut} \
            --mode {wildcards.mmode} \
            --output_1 {output.all_mut} \
            --output_2 {output.min_mut} \
            --log {log}"""


rule count_matrix_split_by_vaf:
    benchmark:
        "workflow/benchmarks/count_matrix_split_by_vaf_{mmode}_{cohort}.tsv"
    input:
        mut=lambda w:config["data"]["mut"][w.cohort],
        mut_sum=lambda w:config["data"]["mut_sum"][w.cohort],
        env="../common/logs/setup_conda_signatures.done"
    conda:
        config["setup"]["Signatures"]
    log:
        "workflow/logs/count_matrix_split_by_vaf_{mmode}_{cohort}.log"
    params:
        min_mut = config["count"]["mutation_min_per_tumor"]*3
    output:
        high_vaf="%s/counts_mutations/counts_mutations_{mmode}_min_mut_high_vaf_{cohort}.tsv" % R_FOLDER,
        low_vaf="%s/counts_mutations/counts_mutations_{mmode}_min_mut_low_vaf_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    threads: 1
    shell:
        """Rscript workflow/scripts/01_count_matrix.R \
            --cohort {wildcards.cohort} \
            --mut {input.mut} \
            --mut_sum {input.mut_sum} \
            --min_mut {params.min_mut} \
            --mode {wildcards.mmode} \
            --split_by_vaf \
            --output_1 {output.high_vaf} \
            --output_2 {output.low_vaf} \
            --log {log}"""
