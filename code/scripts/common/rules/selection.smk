rule selection_samples:
    log:
        "%s/selection_samples.log" % L_FOLDER
    input:
        ["%s/%s" % (D_FOLDER, FILEPATHS[cohort]["bio"]["in_design"]) for cohort in config["selection"]["cohorts"]],
        ["%s/%s" % (D_FOLDER, FILEPATHS[cohort]["cln"]["in_design"]) for cohort in config["selection"]["cohorts"]],
        "../common/logs/setup_conda.done"
    conda:
        "../envs/MetaPrism.yaml"
    output:
        cnt="%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        sam=expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER, cohort=config["selection"]["cohorts"])
    params:
        cohorts = config["selection"]["cohorts"],
        config = "config/config.yaml"
    resources:
        partition = "cpu_short",
        mem_mb = 4000,
        time = "00:15:00"
    threads: 1
    shell:
        """
        Rscript ../common/scripts/get_table_selection.R \
            --cohorts {params.cohorts} \
            --config {params.config} \
            --section selection \
            --output_cnt {output.cnt} \
            --output_sam {output.sam} \
            --log {log}
        """
