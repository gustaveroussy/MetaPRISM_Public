rule draw_heatmap_all_alterations:
    log:
        "%s/draw_heatmap_all_alterations_{direction}.log" % L_FOLDER
    input:
        alterations = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap_all/tables_heatmap_{direction}.xlsx" % R_FOLDER,
        plot = "%s/heatmap_all/heatmap_{direction}.svg" % R_FOLDER,
        plot_paper = "%s/F5_{direction}.svg" % F_FOLDER,
    params:
        width = lambda w: config["heatmap_all"][w.direction]["width"],
        height = lambda w: config["heatmap_all"][w.direction]["height"]
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/02.1_draw_heatmap_all_alterations.R \
            --alterations {input.alterations} \
            --samples {input.samples} \
            --counts {input.counts} \
            --direction {wildcards.direction} \
            --n_cores {threads} \
            --output_tables {output.tables} \
            --output_plot_width {params.width} \
            --output_plot_height {params.height} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --log {log}
        """


rule draw_heatmap_dna:
    log:
        "%s/draw_heatmap_dna_{gene_type}.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        mutations = ["%s/alterations/aggregated_alterations_%s.tsv" % (R_FOLDER, cohort) \
            for cohort in config["data"]["cohorts"]],
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap_dna/tables_{gene_type}.xlsx" % R_FOLDER,
        plot = "%s/heatmap_dna/heatmap_{gene_type}.svg" % R_FOLDER,
        plot_paper = "%s/F3a_{gene_type}.svg" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts_evt = lambda w: [config["heatmap_dna"]["min_counts_evt"][w.gene_type][c] for c in config["data"]["cohorts"]],
        plot_width_one = config["heatmap_dna"]["width_one"],
        plot_height_one = config["heatmap_dna"]["height_one"],
        sample_type = "DNA_T__DNA_N",
        selection_name = "heatmap_dna"
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 4
    shell:
        """
        Rscript workflow/scripts/02.2_draw_heatmap_dna_alterations.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --alterations {input.mutations} \
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
