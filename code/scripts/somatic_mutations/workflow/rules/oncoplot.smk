rule oncoplot_inputs:
    log:
        "workflow/logs/oncoplot_inputs_{cohort}_{tumor_type}.log"
    input:
        mutations = "%s/mutpanning/{cohort}/{tumor_type}_inputs/mutations_{tumor_type}.maf" % R_FOLDER,
        mutpanning_run = "%s/mutpanning/{cohort}/{tumor_type}/SignificanceFiltered/Significance{tumor_type}.txt" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        gen = temp("%s/oncoplots/genes_{cohort}_{tumor_type}.tsv" % R_FOLDER),
        cln = temp("%s/oncoplots/clinical_{cohort}_{tumor_type}.tsv" % R_FOLDER),
        maf = temp("%s/oncoplots/mutations_{cohort}_{tumor_type}.maf" % R_FOLDER),
        lbar = temp("%s/oncoplots/left_bar_data_{cohort}_{tumor_type}.tsv" % R_FOLDER),
        rbar = temp("%s/oncoplots/right_bar_data_{cohort}_{tumor_type}.tsv" % R_FOLDER)
    params:
        vc_selection_list = config["oncoplot"]["vc_selection_list"],
        genes_selection_list = config["oncoplot"]["genes_selection_list"],
        mutpanning_pairs = config["mutpanning"]["pairs"]
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.1_oncoplot_inputs.R \
            --mutations {input.mutations} \
            --cohort {wildcards.cohort} \
            --tumor_type {wildcards.tumor_type} \
            --vc_selection_list {params.vc_selection_list} \
            --genes_selection_list {params.genes_selection_list} \
            --genes_selection_mutpanning_run {input.mutpanning_run} \
            --genes_selection_mutpanning_pairs {params.mutpanning_pairs} \
            --output_gen {output.gen} \
            --output_cln {output.cln} \
            --output_maf {output.maf} \
            --output_lbar {output.lbar} \
            --output_rbar {output.rbar} \
            --log {log}
        """


rule oncoplot_maftools:
    log:
        "workflow/logs/oncoplot_maftools_{cohort}_{tumor_type}.log"
    input:
        gen = "%s/oncoplots/genes_{cohort}_{tumor_type}.tsv" % R_FOLDER,
        cln = "%s/oncoplots/clinical_{cohort}_{tumor_type}.tsv" % R_FOLDER,
        maf = "%s/oncoplots/mutations_{cohort}_{tumor_type}.maf" % R_FOLDER,
        lbar = "%s/oncoplots/left_bar_data_{cohort}_{tumor_type}.tsv" % R_FOLDER,
        rbar = "%s/oncoplots/right_bar_data_{cohort}_{tumor_type}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/oncoplots/oncoplot_{cohort}_{tumor_type}.pdf" % R_FOLDER
    params:
        vc_colors = config["oncoplot"]["vc_colors"]
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.2_oncoplot_maftools.R \
            --input_gen {input.gen} \
            --input_cln {input.cln} \
            --input_maf {input.maf} \
            --input_lbar {input.lbar} \
            --input_rbar {input.rbar} \
            --vc_colors {params.vc_colors} \
            --output {output} \
            --log {log}
        """
