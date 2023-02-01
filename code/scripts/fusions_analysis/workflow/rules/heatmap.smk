rule draw_heatmap:
    log:
        "%s/draw_heatmap_{gene_type}.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=COHORTS_NOT_VAL),
        fusions = ["%s/%s/rna/fusions/%s_annotated_filtered.tsv.gz" % (D_FOLDER, cohort, cohort) \
            for cohort in COHORTS_NOT_VAL],
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap/tables_{gene_type}.xlsx" % R_FOLDER,
        plot = "%s/heatmap/heatmap_{gene_type}.svg" % R_FOLDER,
        plot_paper = "%s/F4b_{gene_type}.svg" % F_FOLDER
    params:
        cohorts = COHORTS_NOT_VAL,
        min_counts_evt = lambda w: [config["heatmap"]["min_counts_evt"][w.gene_type][c] for c in COHORTS_NOT_VAL],
        plot_width_one = config["heatmap"]["width_one"],
        plot_height_one = config["heatmap"]["height_one"],
        sample_type = "RNA_T",
        selection_name = "burden"
    threads: 4
    shell:
        """
        Rscript ../combined_alterations/workflow/scripts/02.2_draw_heatmap_dna_alterations.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --alterations {input.fusions} \
            --counts {input.counts} \
            --min_counts_evt {params.min_counts_evt} \
            --gene_type {wildcards.gene_type} \
            --sample_type {params.sample_type} \
            --selection_name {params.selection_name} \
            --n_cores {threads} \
            --output_tables {output.tables} \
            --output_plot_width_one {params.plot_width_one} \
            --output_plot_height_one {params.plot_height_one} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --log {log}
        """
