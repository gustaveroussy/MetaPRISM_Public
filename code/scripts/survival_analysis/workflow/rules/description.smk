rule desc_dist_survival:
    log:
        "%s/desc_dist_survival.log" % L_FOLDER
    input:
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = config["description"]["cohort"],
        min_count = config["description"]["min_count"],
        width = 8,
        height = 10
    output:
        ["%s/description/dist_%s.pdf" % (R_FOLDER, x) for x in ["survival_times", "counts_drug", "counts_mets"]]
    shell:
        """
        Rscript workflow/scripts/01.1_desc_dist_survival.R \
            --cohort {params.cohort} \
            --min_count {params.min_count} \
            --width {params.width} \
            --height {params.height} \
            --outputs {output} \
            --log {log}
        """
