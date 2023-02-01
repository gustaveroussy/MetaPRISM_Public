rule draw_heatmap:
    log:
        "%s/draw_heatmap.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        mutations = expand("%s/{cohort}/wes/germline_maf/germline_calls_pathogenic.maf.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap/tables_heatmap.xlsx" % R_FOLDER,
        plot = "%s/heatmap/heatmap.svg" % R_FOLDER,
        plot_paper = "%s/F3a_germline.svg" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts_evt = lambda w: [config["heatmap"]["min_counts_evt"][c] for c in config["data"]["cohorts"]],
        plot_width_one = config["heatmap"]["width_one"],
        plot_height_one = config["heatmap"]["height_one"],
        sample_type = "DNA_N",
        selection_name = "plotting"
    threads: 4
    shell:
        """
        Rscript ../combined_alterations/workflow/scripts/02.2_draw_heatmap_dna_alterations.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --alterations {input.mutations} \
            --counts {input.counts} \
            --min_counts_evt {params.min_counts_evt} \
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


rule draw_bar_and_pie_charts:
    log:
        "%s/draw_bar_and_pie_charts.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        mutations = expand("%s/{cohort}/wes/germline_maf/germline_calls_pathogenic.maf.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        pathways = "resources/germline_pathways.txt",
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plots = ["%s/bars_pies/%s.pdf" % (R_FOLDER, x) for x in \
            ["bars_fraction_all", "donut_prism_all","donut_prism_all_dna_repair", "donut_prism_per_tumor_type"]],
        plots_paper = ["%s/%s.svg" % (F_FOLDER, x) for x in ["F3c", "F3d_top", "F3d_bot", "F3e"]]
    params:
        cohorts = config["data"]["cohorts"],
        widths = [300, 250, 125, 500],
        heights = [300, 250, 125, 300]
    threads: 4
    shell:
        """
        python workflow/scripts/01.2_bar_and_pie_charts.py \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --pathways {input.pathways} \
            --widths {params.widths} \
            --heights {params.heights} \
            --outputs {output.plots} \
            --outputs_paper {output.plots_paper} &> {log}
        """
