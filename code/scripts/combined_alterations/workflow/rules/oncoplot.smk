rule oncoplot_inputs:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]]),
        sample_type = "|".join([re.escape(x) for x in ["DNA_T", "RNA_T"]])
    log:
        "%s/oncoplot_inputs_{cohort}_{tumor_type}_{sample_type}.log" % L_FOLDER
    input:
        mutations = "%s/alterations/aggregated_alterations_{cohort}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        # gen = temp("%s/oncoplots/genes_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER),
        # cln = temp("%s/oncoplots/clinical_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER),
        # maf = temp("%s/oncoplots/mutations_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER),
        # lbar = temp("%s/oncoplots/left_bar_data_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER)
        gen = "%s/oncoplots/genes_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        cln = "%s/oncoplots/clinical_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        maf = "%s/oncoplots/mutations_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        lbar = "%s/oncoplots/left_bar_data_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER
    params:
        vc_selection_list = lambda w: config["oncoplot"]["vc_selection_list"][w.sample_type],
        genes_selection_mode = "min_freq",
        genes_selection_list = config["oncoplot"]["genes_selection_list"],
        genes_selection_min_freq = lambda w: config["oncoplot"]["genes_selection_min_freq"][w.cohort],
        algorithms = lambda w: config["oncoplot"]["algorithms"][w.sample_type]
    threads: 1
    shell:
        """
        Rscript ../somatic_mutations/workflow/scripts/02.1_oncoplot_inputs.R \
            --mutations {input.mutations} \
            --cohort {wildcards.cohort} \
            --tumor_type {wildcards.tumor_type} \
            --sample_type {wildcards.sample_type} \
            --algorithms {params.algorithms} \
            --vc_selection_list {params.vc_selection_list} \
            --genes_selection_mode {params.genes_selection_mode} \
            --genes_selection_list {params.genes_selection_list} \
            --genes_selection_min_freq {params.genes_selection_min_freq} \
            --output_gen {output.gen} \
            --output_cln {output.cln} \
            --output_maf {output.maf} \
            --output_lbar {output.lbar} \
            --log {log}
        """


rule oncoplot_maftools:
    log:
        "workflow/logs/oncoplot_maftools_{cohort}_{tumor_type}_{sample_type}.log"
    input:
        gen = "%s/oncoplots/genes_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        cln = "%s/oncoplots/clinical_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        maf = "%s/oncoplots/mutations_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER,
        lbar = "%s/oncoplots/left_bar_data_{cohort}_{tumor_type}_{sample_type}.tsv" % R_FOLDER
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/oncoplots/oncoplot_{cohort}_{tumor_type}_{sample_type}.pdf" % R_FOLDER
    params:
        show_samples = lambda w: "no" if w.cohort=="tcga" else "yes",
        vc_colors = config["oncoplot"]["vc_colors"],
        font_size = 0.4
    threads: 1
    shell:
        """
        Rscript ../somatic_mutations/workflow/scripts/02.2_oncoplot_maftools.R \
            --input_gen {input.gen} \
            --input_cln {input.cln} \
            --input_maf {input.maf} \
            --input_lbar {input.lbar} \
            --show_samples {params.show_samples} \
            --vc_colors {params.vc_colors} \
            --font_size {params.font_size} \
            --output {output} \
            --log {log}
        """
