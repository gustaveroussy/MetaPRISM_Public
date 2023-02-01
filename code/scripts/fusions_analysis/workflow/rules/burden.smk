rule burden_tables_sel:
    log:
        "workflow/logs/burden_tables_sel.log"
    input:
        fusions = expand("%s/tables/{cohort}/fusions_preprocessed.tsv" % R_FOLDER, cohort=COHORTS),
        samples = expand("%s/tables/{cohort}/selection_samples.tsv" % R_FOLDER, cohort=COHORTS),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    params:
        cohorts = COHORTS
    output:
        expand("%s/tables/{cohort}/burden_table_sel.tsv" % R_FOLDER, cohort=COHORTS)
    shell:
        """Rscript workflow/scripts/01.1_burden_tables.R \
            --input_cohorts {params.cohorts} \
            --input_samples {input.samples} \
            --input_fusions {input.fusions} \
            --input_counts {input.counts} \
            --output {output} \
            --log {log}"""

rule plot_burden_sel:
    log:
        "workflow/logs/plot_burden_sel.log"
    input:
        tables = rules.burden_tables_sel.output,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    params:
        cohorts = COHORTS,
        output_width = config["burden"]["sel"]["output_width"],
        output_height = config["burden"]["sel"]["output_height"]
    output:
        "%s/plots/plot_burden_tumor_types_sel.pdf" % R_FOLDER
    shell:
        """Rscript workflow/scripts/01.2_plot_burden.R \
            --input_cohorts {params.cohorts} \
            --input_tables {input.tables} \
            --output_width {params.output_width} \
            --output_height {params.output_height} \
            --output {output} \
            --log {log}"""


rule burden_tables_all:
    log:
        "workflow/logs/burden_tables_all.log"
    input:
        fusions = expand("%s/tables/{cohort}/fusions_preprocessed.tsv" % R_FOLDER, cohort=COHORTS),
        samples = expand("%s/tables/{cohort}/selection_samples.tsv" % R_FOLDER, cohort=COHORTS),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    params:
        cohorts = COHORTS
    output:
        expand("%s/tables/{cohort}/burden_table_all.tsv" % R_FOLDER, cohort=COHORTS)
    shell:
        """Rscript workflow/scripts/01.1_burden_tables.R \
            --input_cohorts {params.cohorts} \
            --input_samples {input.samples} \
            --input_fusions {input.fusions} \
            --input_counts {input.counts} \
            --output {output} \
            --log {log}"""

rule plot_burden_all:
    log:
        "workflow/logs/plot_burden_all.log"
    input:
        tables = rules.burden_tables_all.output,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    params:
        cohorts = COHORTS,
        output_width = config["burden"]["sel"]["output_width"],
        output_height = config["burden"]["sel"]["output_height"]
    output:
        "%s/plots/plot_burden_tumor_types_all.pdf" % R_FOLDER
    shell:
        """Rscript workflow/scripts/01.2_plot_burden.R \
            --input_cohorts {params.cohorts} \
            --input_tables {input.tables} \
            --output_width {params.output_width} \
            --output_height {params.output_height} \
            --output {output} \
            --log {log}"""
