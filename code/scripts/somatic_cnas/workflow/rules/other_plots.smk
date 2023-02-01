rule draw_cna_summary:
    log:
        "%s/draw_cna_summary.log" % L_FOLDER
    input:
        sam_tables = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        cln_tables = expand("%s/{cohort}/clinical/curated/cln_{cohort}_in_design_curated.tsv" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        cna_tables = expand("%s/{cohort}/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plots = ["%s/other_plots/density_ploidies.pdf" % R_FOLDER,
                 "%s/other_plots/boxplot_losses_and_gains_all.pdf",
                 "%s/other_plots/boxplot_llg_mlg_and_loh_cnloh_prism_wgd.pdf",
                 "%s/other_plots/boxplot_hlg_and_hd_prism_wgd.pdf",
                 "%s/other_plots/boxplot_chr_gains_and_chr_losses_prism_wgd.pdf",
                 "%s/other_plots/boxplot_ogs_and_tsgs_prism_wgd.pdf"],
        plots_paper = ["%s/%s" % (F_FOLDER, x) for x \
            in ["F2c.eps", "F2d.eps", "FS11a.pdf", "FS11b.pdf", "FS11c.pdf", "FS11d.pdf"]]
    params:
        cohorts = config["data"]["cohorts"],
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.2_draw_cna_summary.R \
            --cohorts {params.cohorts} \
            --sam_tables {input.sam_tables} \
            --cna_tables {input.cna_tables} \
            --cln_tables {input.cln_tables} \
            --output_plots {output.plots} \
            --output_plots_paper {output.plots_paper} \
            --log {log}
        """


rule draw_cna_barplots:
    log:
        "%s/draw_cna_barplots.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        cln_tables = expand("%s/{cohort}/clinical/curated/cln_{cohort}_in_design_curated.tsv" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        cna_tables = expand("%s/{cohort}/wes/somatic_cna/somatic_calls.tsv.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plot = "%s/other_plots/barplots_main_genes.pdf" % R_FOLDER,
        plot_paper = "%s/FS11d-e.pdf" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.3_draw_cna_barplots.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --cln_tables {input.cln_tables} \
            --cna_tables {input.cna_tables} \
            --counts {input.counts} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --log {log}
        """
