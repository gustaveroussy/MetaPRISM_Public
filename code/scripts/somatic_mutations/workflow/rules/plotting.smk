rule burden_prism_tumor_types:
    log:
        "%s/burden_prism_tumor_types.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER, cohort=config["data"]["cohorts"]),
        mutations = [config["data"]["mut"]["sum"][cohort] for cohort in config["data"]["cohorts"]],
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        target_bed = config["data"]["target_bed"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plot="%s/burden/cumulative_scatter_plot_prism_tumor_types.pdf" % R_FOLDER,
        plot_paper="%s/F2a.svg" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts = 10,
        margin_top = 70,
        plot_width = 580,
        plot_height = 250
    threads: 1
    shell:
        """
        Rscript workflow/scripts/03.1_plot_burden.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --target_bed {input.target_bed} \
            --counts {input.counts} \
            --min_counts {params.min_counts} \
            --margin_top {params.margin_top} \
            --output_width {params.plot_width} \
            --output_height {params.plot_height} \
            --output {output.plot} \
            --output_paper {output.plot_paper} \
            --log {log}
        """


rule burden_all_tumor_types:
    log:
        "%s/burden_all_tumor_types.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER, cohort=config["data"]["cohorts"]),
        mutations = [config["data"]["mut"]["sum"][cohort] for cohort in config["data"]["cohorts"]],
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        target_bed = config["data"]["target_bed"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plot="%s/burden/cumulative_scatter_plot_all_tumor_types.pdf" % R_FOLDER,
        plot_paper="%s/F2a_1.pdf" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts = 10,
        margin_top = 120,
        plot_width = 1600,
        plot_height = 300
    threads: 1
    shell:
        """
        Rscript workflow/scripts/03.1_plot_burden.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --target_bed {input.target_bed} \
            --counts {input.counts} \
            --show_pvals \
            --min_counts {params.min_counts} \
            --margin_top {params.margin_top} \
            --use_all_samples \
            --output_width {params.plot_width} \
            --output_height {params.plot_height} \
            --output {output.plot} \
            --output_paper {output.plot_paper} \
            --log {log}
        """


rule barplots_recurrence_mutations:
    log:
        "%s/barplots_recurrence_{selection_mut}_mutations_cby_{color_by}_gby_{group_by}_in_{tumor_type}_of_{cohort}.log" % L_FOLDER
    input:
        mutations_all = lambda w: config["data"]["mut"]["all"][w.cohort],
        mutations_ann = lambda w: config["data"]["mut"]["ann"][w.cohort],
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        recurrence_threshold = lambda w: config["recurrence"]["threshold"][w.group_by][w.selection_mut][w.cohort]
    output:
        "%s/recurrence/barplots_{selection_mut}_mutations_cby_{color_by}_gby_{group_by}_in_{tumor_type}_of_{cohort}.pdf" % R_FOLDER
    shell:
        """
        python workflow/scripts/05.1_barplots_recurrence.py \
            --cohort {wildcards.cohort} \
            --mutations_all {input.mutations_all} \
            --mutations_ann {input.mutations_ann} \
            --selection_mut {wildcards.selection_mut} \
            --tumor_type {wildcards.tumor_type} \
            --color_by {wildcards.color_by} \
            --group_by {wildcards.group_by} \
            --recurrence_threshold {params.recurrence_threshold} \
            --output {output} &> {log}
        """
