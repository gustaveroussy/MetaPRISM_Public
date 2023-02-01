rule draw_heatmap_chr_arm:
    log:
        "%s/draw_heatmap_chr_arm_{alterations}.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        cna_chr = expand("%s/{cohort}/wes/somatic_cna/somatic_calls_per_chr_arm.tsv" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        cna_wgd = expand("%s/{cohort}/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap_chr_arm/tables_{alterations}.xlsx" % R_FOLDER,
        plot = "%s/heatmap_chr_arm/heatmap_{alterations}.pdf" % R_FOLDER,
        plot_paper = "%s/FS10_{alterations}.pdf" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts_evt = lambda w: [config["heatmap_chr_arm"]["min_counts_evt"][w.alterations][c] for c in config["data"]["cohorts"]],
        plot_width_one = config["heatmap_chr_arm"]["width_one"],
        plot_height_one = config["heatmap_chr_arm"]["height_one"]
    threads: 4
    shell:
        """
        Rscript workflow/scripts/02.1_draw_heatmap_chr_arm.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --cna_chr {input.cna_chr} \
            --cna_wgd {input.cna_wgd} \
            --counts {input.counts} \
            --min_counts_evt {params.min_counts_evt} \
            --alterations {wildcards.alterations} \
            --n_cores {threads} \
            --output_tables {output.tables} \
            --output_plot_width_one {params.plot_width_one} \
            --output_plot_height_one {params.plot_height_one} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --log {log}
        """
