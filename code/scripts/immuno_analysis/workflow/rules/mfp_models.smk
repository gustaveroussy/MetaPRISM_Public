rule mfp_tcga_signatures_mfp_model:
    input:
        signatures = "%s/tcga/tables/signatures_bagaev_2021_mfp.tsv" % R_FOLDER,
        annotation = config["resources"]["bagaev_2021"]["annotation"],
        bio = "%s/tcga/clinical/curated_other/bio_tcga_all_curated.tsv" % D_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_tcga_signatures_mfp_model.log" % L_FOLDER
    output:
        signatures = "%s/tcga/tables/signatures_mfp_model.tsv" % R_FOLDER,
        annotation = "%s/tcga/tables/annotation_mfp_model.tsv" % R_FOLDER
    threads: 1
    shell:
        """
        python workflow/scripts/03.2_mfp_tcga_signatures_mfp_model.py \
            --input_signatures {input.signatures} \
            --input_annotation {input.annotation} \
            --input_bio {input.bio} \
            --output_signatures {output.signatures} \
            --output_annotation {output.annotation} &> {log}
        """


rule mfp_model:
    benchmark:
        "%s/mfp_model_{cohort}_{model_name}.tsv" % B_FOLDER
    input:
        signatures = get_mfp_model_signatures,
        annotation = get_mfp_model_annotation,
        models = get_input_models_mfp_model,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_model_{cohort}_{model_name}.log" % L_FOLDER
    params:
        grids_yaml = config["mfp_models"]["grids_yaml"],
        model_name = lambda w: w.model_name,
        random_state = config["random_state"]
    output:
        "%s/{cohort}/mfp_models/{model_name}/pipeline_train.joblib" % R_FOLDER
    threads: get_threads_mfp_model
    resources:
        mem_mb = get_mem_mb_mfp_model,
        partition = get_partition_mfp_model,
        time = get_time_mfp_model
    shell:
        """
        python workflow/scripts/03.4_mfp_models_train.py \
            --input_signatures {input.signatures} \
            --input_annotation {input.annotation} \
            --grids_yaml {params.grids_yaml} \
            --model_name {params.model_name} \
            --random_state {params.random_state} \
            --threads {threads} \
            --output {output} &> {log}
        """


rule mfp_predict:
    benchmark:
        "%s/mfp_predict_{cohort}_{model_name}.tsv" % B_FOLDER
    input:
        samples = get_mfp_predict_samples,
        signatures = get_mfp_predict_signatures,
        model = get_mfp_predict_model,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_predict_{cohort}_{model_name}.log" % L_FOLDER
    output:
        "%s/{cohort}/tables/mfp_subtypes_predicted_{model_name}.tsv" % R_FOLDER
    params:
        input_annotation = get_mfp_predict_annotation
    threads: 1
    resources:
        mem_mb = 16000,
        partition = "cpu_short",
        time = "01:00:00"
    shell:
        """
        python workflow/scripts/03.5_mfp_models_predict.py \
            --input_samples {input.samples} \
            --input_signatures {input.signatures} \
            --input_annotation {params.input_annotation} \
            --input_model {input.model} \
            --output {output} &> {log}
        """


rule mfp_extract:
    input:
        samples = get_mfp_predict_samples,
        signatures = get_mfp_predict_signatures,
        model = get_mfp_extract_model,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_extract_{cohort}.log" % L_FOLDER
    output:
        "%s/{cohort}/tables/signatures_mfp_model_preprocessed.tsv" % R_FOLDER
    params:
        input_annotation = get_mfp_predict_annotation
    threads: 1
    resources:
        mem_mb = 16000,
        partition = "cpu_short",
        time = "01:00:00"
    shell:
        """
        python workflow/scripts/03.6_mfp_models_extract_preprocessed_signatures.py \
            --input_samples {input.samples} \
            --input_signatures {input.signatures} \
            --input_annotation {params.input_annotation} \
            --input_model {input.model} \
            --output {output} &> {log}
        """


rule mfp_heatmap:
    input:
        signatures = "%s/{cohort}/tables/signatures_mfp_model_preprocessed.tsv" % R_FOLDER,
        cln = "%s/{cohort}/clinical/curated/cln_{cohort}_in_design_curated.tsv" % D_FOLDER,
        bio = "%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER,
        subtypes = "%s/{cohort}/tables/mfp_subtypes_predicted_LogisticRegression.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_heatmap_{cohort}.log" % L_FOLDER
    output:
        directory("%s/{cohort}/plots" % R_FOLDER)
    threads: 1
    resources:
        mem_mb = 4000,
        partition = "cpu_short",
        time = "01:00:00"
    shell:
        """
        Rscript workflow/scripts/03.7_mfp_heatmap.R \
            --input_signatures {input.signatures} \
            --input_cln {input.cln} \
            --input_bio {input.bio} \
            --input_subtypes {input.subtypes} \
            --output {output} \
            --log {log}
        """

rule mfp_heatmap_summary:
    log:
        "%s/mfp_heatmap_summary.log" % L_FOLDER
    input:
        samples = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        subtypes = expand("%s/{cohort}/tables/mfp_subtypes_predicted_LogisticRegression.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        tables = "%s/heatmap/tables.xlsx" % R_FOLDER,
        plot = "%s/heatmap/heatmap.svg" % R_FOLDER,
        plot_paper = "%s/F4a.svg" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"],
        min_counts_evt = [1 for _ in config["data"]["cohorts"]],
        plot_width_one = config["heatmap"]["width_one"],
        plot_height_one = config["heatmap"]["height_one"],
        sample_type = "RNA_T",
        selection_name = "mfp_model",
        gene_type = "tmes"
    threads: 4
    shell:
        """
        Rscript ../combined_alterations/workflow/scripts/02.2_draw_heatmap_dna_alterations.R \
            --cohorts {params.cohorts} \
            --samples {input.samples} \
            --alterations {input.subtypes} \
            --counts {input.counts} \
            --min_counts_evt {params.min_counts_evt} \
            --gene_type {params.gene_type} \
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
