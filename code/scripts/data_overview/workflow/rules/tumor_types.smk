rule tt_donut_plots:
    log:
        "%s/tt_donut_plots.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohorts=config["tumor_types"]["cohorts"],
        width=400,
        height=400
    output:
        out=expand("%s/tumor_types/donut_plot_tumor_types_{cohort}.svg" % R_FOLDER,
            cohort=config["tumor_types"]["cohorts"]),
        out_paper=expand("%s/F{x}" % F_FOLDER, x=["1a_left.svg","S3c.pdf","S3a.pdf"])
    shell:
        """
        python -u workflow/scripts/01.1_tt_donut_plots.py \
            --cohorts {params.cohorts} \
            --width {params.width} \
            --height {params.height} \
            --outputs {output.out} \
            --outputs_paper {output.out_paper} &> {log}
        """

rule tt_barplots:
    log:
        "%s/tt_barplots.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohorts=config["tumor_types"]["cohorts"],
        min_count=config["tumor_types"]["min_count"],
        tt_drop=config["tumor_types"]["tt_drop"],
        width=1000,
        height=500
    output:
        samples=expand("%s/tumor_types/bar_chart_tumor_types_sample_types_{cohort}.svg" % R_FOLDER,
            cohort=config["tumor_types"]["cohorts"]),
        all="%s/tumor_types/bar_chart_tumor_types_sample_types_all_cohorts.svg" % R_FOLDER,
        rare=expand("%s/tumor_types/bar_chart_rare_histological_types_{cohort}.pdf" % R_FOLDER,
            cohort=["prism", "met500"]),
        paper=expand("%s/F{x}" % F_FOLDER, x=["S1.pdf","1b.svg","S3d.pdf", "S3b.pdf"])
    shell:
        """
        python workflow/scripts/01.2_tt_barplots.py \
            --cohorts {params.cohorts} \
            --min_count {params.min_count} \
            --tt_drop {params.tt_drop} \
            --width {params.width} \
            --height {params.height} \
            --outputs_samples {output.samples} \
            --output_all {output.all} \
            --outputs_rare {output.rare} \
            --outputs_paper {output.paper} &> {log}
        """
