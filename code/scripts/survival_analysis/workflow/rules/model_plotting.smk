rule plot_models_coeffs:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
        model = "|".join([re.escape(x) for x in MODELS]),
    benchmark:
        "%s/plot_models_coeffs_{cohort}_{samples}_{features}_{selection}_{model}.tsv" % B_FOLDER
    log:
        "%s/plot_models_coeffs_{cohort}_{samples}_{features}_{selection}_{model}.log" % L_FOLDER
    input:
        met="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.pooled_ax_rep.tsv.gz" % R_FOLDER,
        cov="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.pooled_ax_rep.tsv.gz" % R_FOLDER,
        cov_final="%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER,
        cov_remov="%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.removed.tsv.gz" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        threshold = 0.2,
        ci_conf = 0.95,
        outcome = "survival"
    threads: 1
    output:
        forester="%s/plots_{cohort}/{samples}/{features}/forester_coeffs_{selection}_{model}.pdf" % R_FOLDER,
        pieplots="%s/plots_{cohort}/{samples}/{features}/pieplots_coeffs_{selection}_{model}.pdf" % R_FOLDER
    resources:
        partition = "cpu_med",
        time = "00:10:00"
    shell:
        """
        Rscript ../common/scripts/model_plots_coeffs.R \
            --input_met {input.met} \
            --input_cov {input.cov} \
            --cov_final {input.cov_final} \
            --cov_remov {input.cov_remov} \
            --threshold {params.threshold} \
            --ci_conf {params.ci_conf} \
            --outcome {params.outcome} \
            --output_forester {output.forester} \
            --output_pieplots {output.pieplots} \
            --log {log}
        """


rule plot_metrics_per_selection:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
    benchmark:
        "%s/plot_metrics_{cohort}_{samples}_{time_horizon}_{selection}_results.tsv" % B_FOLDER
    log:
        "%s/plot_metrics_{cohort}_{samples}_{time_horizon}_{selection}_results.log" % L_FOLDER
    input:
        get_input_plot_metrics,
        "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        names_features = FEATURES,
        names_models = MODELS,
        names_selections = lambda w: [w.selection],
        group_by = "Model",
        dir_data = lambda w: "%s/data_{cohort}/sub_features/%s" % (R_FOLDER, w.samples),
        dir_models = lambda w: "%s/models_{cohort}/%s" % (R_FOLDER, w.samples),
        drop_warn = "true",
        drop_error = "true",
        width_met_one = 140,
        height_met = 1000
    threads: 1
    output:
        met_bscore="%s/plots_{cohort}/{samples}/boxplot_bscore_{time_horizon}_{selection}.pdf" % R_FOLDER,
        met_cindex="%s/plots_{cohort}/{samples}/boxplot_cindex_{time_horizon}_{selection}.pdf" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 10000,
        time = "00:20:00"
    shell:
        """
        python -u workflow/scripts/04.1_plot_metrics.py \
            --names_features {params.names_features} \
            --names_models {params.names_models} \
            --names_selections {params.names_selections} \
            --group_by {params.group_by} \
            --dir_data {params.dir_data} \
            --dir_models {params.dir_models} \
            --time_horizon {wildcards.time_horizon} \
            --drop_warn {params.drop_warn} \
            --drop_error {params.drop_error} \
            --width_met_one {params.width_met_one} \
            --height_met {params.height_met} \
            --output_met_bscore {output.met_bscore} \
            --output_met_cindex {output.met_cindex} &> {log}
        """


rule plot_discretized_risk_score:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        cohort_pred = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
        model = "|".join([re.escape(x) for x in MODELS]),
    benchmark:
        "%s/plot_discretize_risk_score_{cohort}_{samples}_{features}_{selection}_{model}_{cohort_pred}.tsv" % B_FOLDER
    log:
        "%s/plot_discretize_risk_score_{cohort}_{samples}_{features}_{selection}_{model}_{cohort_pred}.log" % L_FOLDER
    input:
        dats = get_input_dats_plot_discretized_risk_score,
        cov_model = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER,
        cov_dats = "%s/data_{cohort_pred}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER,
        model = "%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.pooled_ax_rep.tsv.gz" % R_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        name_model = get_name_model_plot_discretized_risk_score
    threads: 1
    output:
        densities = "%s/plots_{cohort}/{samples}/{features}/densities_6_month_{selection}_{model}_{cohort_pred}.pdf" % R_FOLDER,
        km_curves = "%s/plots_{cohort}/{samples}/{features}/km_curves_risk_{selection}_{model}_{cohort_pred}.pdf" % R_FOLDER,
        risks = "%s/models_{cohort}/{samples}/{features}/{selection}_{model}/predicted_risk_scores_{cohort_pred}.tsv.gz" % R_FOLDER,
    resources:
        partition = "cpu_med",
        mem_mb = 4000,
        time = "00:20:00"
    shell:
        """
        Rscript workflow/scripts/04.3_plot_discretized_risk_score.R \
            --name_cohort {wildcards.cohort_pred} \
            --name_samples {wildcards.samples} \
            --input_dats {input.dats} \
            --input_model {input.model} \
            --input_cov_dats {input.cov_dats} \
            --input_cov_model {input.cov_model} \
            --output_densities {output.densities} \
            --output_km_curves {output.km_curves} \
            --output_risks {output.risks} \
            --name_model {params.name_model} \
            --log {log}
        """
