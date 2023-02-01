rule bootstrap_cross_validation_indices:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES])
    benchmark:
        "%s/bootstrap_cross_validation_indices_{cohort}_{samples}_{features}.log" % B_FOLDER
    log:
        "%s/bootstrap_cross_validation_indices_{cohort}_{samples}_{features}.log" % L_FOLDER
    input:
        cov = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        dat_folder = "%s/data_{cohort}/sub_features/{samples}/{features}/processed" % R_FOLDER,
        n_boot = config["models"]["fitting"]["n_boot"],
        n_xval = config["models"]["fitting"]["n_xval"],
        k_cv = config["models"]["fitting"]["k_cv"],
        seed = config["models"]["fitting"]["seed"],
    threads: 1
    output:
        boot="%s/models_{cohort}/{samples}/{features}/boot_indices.tsv.gz" % R_FOLDER,
        xval="%s/models_{cohort}/{samples}/{features}/xval_indices.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 4000,
        time = "00:15:00"
    shell:
        """
        Rscript workflow/scripts/03.1_bootstrap_cross_validation_indices.R \
            --dat_folder {params.dat_folder} \
            --n_boot {params.n_boot} \
            --n_xval {params.n_xval} \
            --k_cv {params.k_cv} \
            --seed {params.seed} \
            --output_boot {output.boot} \
            --output_xval {output.xval} \
            --log {log}
        """


rule run_survival_model:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
        model = "|".join([re.escape(x) for x in MODELS]),
    benchmark:
        "%s/run_survival_model_{cohort}_{samples}_{features}_{table}_{selection}_{model}.tsv" % B_FOLDER
    log:
        "%s/run_survival_model_{cohort}_{samples}_{features}_{table}_{selection}_{model}.log" % L_FOLDER
    input:
        dat = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/data.{table}.tsv.gz" % R_FOLDER,
        cov = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER,
        ind = get_input_indices_run_survival,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    threads: config["models"]["fitting"]["threads"]
    # threads: 1
    params:
        repeat_mode = config["models"]["fitting"]["repeat_mode"]
    output:
        # met=temp("%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.{table}.tsv.gz" % R_FOLDER),
        # cov=temp("%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.{table}.tsv.gz" % R_FOLDER)
        met="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.{table}.tsv.gz" % R_FOLDER,
        cov="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.{table}.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 1000*config["models"]["fitting"]["threads"],
        # mem_mb = 4000,
        time = "01:00:00"
    shell:
        """
        Rscript workflow/scripts/03.2_run_survival_model_repeat.R \
            --name_features {wildcards.features} \
            --name_table {wildcards.table} \
            --name_selection {wildcards.selection} \
            --name_model {wildcards.model} \
            --input_dat {input.dat} \
            --input_cov {input.cov} \
            --indices {input.ind} \
            --repeat_mode {params.repeat_mode} \
            --verbose \
            --n_cores {threads} \
            --output_met {output.met} \
            --output_cov {output.cov} \
            --log {log}
        """


rule pool_model_results_ax_imputations:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
        model = "|".join([re.escape(x) for x in MODELS])
    benchmark:
        "%s/pool_model_results_ax_imputations_{cohort}_{samples}_{features}_{selection}_{model}.tsv" % B_FOLDER
    log:
        "%s/pool_model_results_ax_imputations_{cohort}_{samples}_{features}_{selection}_{model}.log" % L_FOLDER
    input:
        get_input_pool_ax_imputations
    conda:
        config["setup"]["MetaPrism"]
    params:
        n_imputations = config["models"]["preprocess"]["n_imputations"],
        run_results = "%s/models_{cohort}/{samples}/{features}/{selection}_{model}" % R_FOLDER
    threads: 1
    output:
        met="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.pooled_ax_imp.tsv.gz" % R_FOLDER,
        cov="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.pooled_ax_imp.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 2000,
        time = "00:10:00"
    shell:
        """
        python ../common/scripts/model_pool_results_ax_imputations.py \
            --run_results {params.run_results} \
            --n_imputations {params.n_imputations} \
            --output_met {output.met} \
            --output_cov {output.cov} &> {log}
        """


rule pool_model_results_ax_repeats:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES]),
        selection = "|".join([re.escape(x) for x in SELECTIONS]),
        model = "|".join([re.escape(x) for x in MODELS])
    benchmark:
        "%s/pool_model_results_ax_repeats_{cohort}_{samples}_{features}_{selection}_{model}.tsv" % B_FOLDER
    log:
        "%s/pool_model_results_ax_repeats_{cohort}_{samples}_{features}_{selection}_{model}.log" % L_FOLDER
    input:
        met="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.pooled_ax_imp.tsv.gz" % R_FOLDER,
        cov="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.pooled_ax_imp.tsv.gz" % R_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        ci_types = ["basic", "percentile"],
        ci_conf = 0.95,
        drop_warn = "true",
        drop_error = "true"
    threads: 1
    output:
        met="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/mets.pooled_ax_rep.tsv.gz" % R_FOLDER,
        cov="%s/models_{cohort}/{samples}/{features}/{selection}_{model}/covs.pooled_ax_rep.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 2000,
        time = "00:10:00"
    shell:
        """
        Rscript ../common/scripts/model_pool_results_ax_repeats.R \
            --input_met {input.met} \
            --input_cov {input.cov} \
            --ci_types {params.ci_types} \
            --ci_conf {params.ci_conf} \
            --drop_warn {params.drop_warn} \
            --drop_error {params.drop_error} \
            --output_met {output.met} \
            --output_cov {output.cov} \
            --log {log}
        """
