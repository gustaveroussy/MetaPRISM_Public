rule upset_plots_validation:
    log:
        "workflow/logs/upset_plots_validation_{cohort}.log"
    input:
        fusions="%s/{cohort}/rna/fusions/{cohort}_annotated.tsv.gz" % D_FOLDER,
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        samples="%s/validation/upset_samples_{cohort}.pdf" % R_FOLDER,
        wo_breakpoints_all="%s/validation/upset_callers_wo_breakpoints_all_{cohort}.pdf" % R_FOLDER,
        wo_breakpoints_com="%s/validation/upset_callers_wo_breakpoints_com_{cohort}.pdf" % R_FOLDER,
        w_breakpoints_all="%s/validation/upset_callers_w_breakpoints_all_{cohort}.pdf" % R_FOLDER,
        w_breakpoints_com="%s/validation/upset_callers_w_breakpoints_com_{cohort}.pdf" % R_FOLDER
    params:
        algos_tcga=config["data"]["aggregate"]["tcga"]["algos"],
        algos_tcga_validation=config["data"]["aggregate"]["tcga_validation"]["algos"]
    resources:
        mem_mb=8000,
        partition="shortq",
        time_min=20
    threads: 1
    shell:
        """Rscript workflow/scripts/00.5_upset_plots_validation.R \
            --fusions {input.fusions} \
            --cohort {wildcards.cohort} \
            --algos_tcga {params.algos_tcga} \
            --algos_tcga_validation {params.algos_tcga_validation} \
            --output_samples {output.samples} \
            --output_wo_breakpoints_all {output.wo_breakpoints_all} \
            --output_wo_breakpoints_com {output.wo_breakpoints_com} \
            --output_w_breakpoints_all {output.w_breakpoints_all} \
            --output_w_breakpoints_com {output.w_breakpoints_com} \
            --log {log}
        """


rule find_best_filtering_criteria:
    log:
        "workflow/logs/find_best_filtering_criteria_use_whitelist_{whitelist}_use_breakpoints_{breakpoints}.log"
    input:
        fusions_tcga="%s/tcga/rna/fusions/tcga_annotated.tsv.gz" % D_FOLDER,
        fusions_tcga_validation="%s/tcga_validation/rna/fusions/tcga_validation_annotated.tsv.gz" % D_FOLDER,
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
       "%s/validation/table_scores_use_whitelist_{whitelist}_use_breakpoints_{breakpoints}.tsv" % R_FOLDER
    params:
        algos_tcga=config["data"]["aggregate"]["tcga"]["algos"],
        algos_tcga_validation=["arriba", "ericscript", "pizzly", "starfusion"]
    resources:
        mem_mb=4000,
        partition="shortq",
        time_min=30
    threads: 16
    shell:
        """python workflow/scripts/00.6_find_best_filtering_criteria.py \
            --fusions_tcga {input.fusions_tcga} \
            --fusions_tcga_validation {input.fusions_tcga_validation} \
            --algos_tcga {params.algos_tcga} \
            --algos_tcga_validation {params.algos_tcga_validation} \
            --use_whitelist {wildcards.whitelist} \
            --use_breakpoints {wildcards.breakpoints} \
            --n_jobs {threads} \
            --output_table {output} &> {log}
        """
