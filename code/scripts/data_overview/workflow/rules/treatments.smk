rule trt_sankeys:
    log:
        "%s/trt_sankeys.log" % L_FOLDER
    input:
        drug_table="%s/%s" % (D_FOLDER, FILEPATHS["resources"]["drug_table"]),
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        widths=[390, 500, 500],
        heights=[1275, 650, 925]
    output:
        "%s/treatments/sankeys/sankey_tumor_type_drug_name_prism.svg" % R_FOLDER,
        "%s/treatments/sankeys/sankey_tumor_type_drug_class_prism.svg" % R_FOLDER,
        "%s/treatments/sankeys/sankey_tumor_type_drug_class_prism_rare.svg" % R_FOLDER,
        "%s/treatments/side_barplot_platinum.pdf" % R_FOLDER,
        "%s/treatments/stacked_barplot_platinum.pdf" % R_FOLDER,
        "%s/F2b_mid.svg" % F_FOLDER,
    shell:
        """
        python workflow/scripts/02.1_trt_sankeys.py \
            --cohort {params.cohort} \
            --drug_table {input.drug_table} \
            --widths {params.widths} \
            --heights {params.heights} \
            --outputs {output} &> {log}
        """


rule trt_barplots_description:
    log:
        "%s/trt_barplots_description_{tumor_type}.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        modes=["All_Therapies", "Targeted_Therapies"]
    output:
        all_therapies="%s/treatments/barplots/{tumor_type}/trt_classes.pdf" % R_FOLDER,
        targeted_therapies="%s/treatments/barplots/{tumor_type}/trt_classes_targeted_therapies.pdf" % R_FOLDER
    shell:
        """
        python workflow/scripts/02.2_trt_barplots_description.py \
            --cohort {params.cohort} \
            --tumor_type "{wildcards.tumor_type}" \
            --modes {params.modes} \
            --outputs {output} &> {log}
        """



rule trt_barplots_resistance:
    log:
        "%s/trt_barplots_resistance_{tumor_type}.log" % L_FOLDER
    input:
        drug_table="%s/%s" % (D_FOLDER, FILEPATHS["resources"]["drug_table"]),
        alts_table=config["data"]["alts_table"],
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        widths=[1400, 1000],
        heights=[500, 650]
    output:
        histologies=["%s/treatments/barplots/{tumor_type}/trt_%s_vs_histologies_prism.pdf" % (R_FOLDER,x)
            for x in ["name", "class"]],
        alterations=["%s/treatments/barplots/{tumor_type}/trt_%s_vs_alterations_prism.pdf" % (R_FOLDER,x)
            for x in ["name", "class"]]
    shell:
        """
        python workflow/scripts/02.2_trt_barplots_resistance.py \
            --cohort {params.cohort} \
            --tumor_type "{wildcards.tumor_type}" \
            --drug_table {input.drug_table} \
            --alts_table {input.alts_table} \
            --widths {params.widths} \
            --heights {params.heights} \
            --output_histologies {output.histologies} \
            --output_alterations {output.alterations} &> {log}
        """


rule trt_forester_1:
    log:
        "%s/trt_forester_1_{drug}_{tumor_type}.log" % L_FOLDER
    input:
        drug_table="%s/%s" % (D_FOLDER, FILEPATHS["resources"]["drug_table"]),
        sigs_table=config["data"]["sigs_table"],
        alts_table=config["data"]["alts_table"],
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        combination="union",
        plot_mode=1,
        drug_classes=lambda w: config["treatments"]["forester_1"][w.drug]["drug_classes"],
        drug_names=lambda w: config["treatments"]["forester_1"][w.drug]["drug_names"],
        sigs_names=lambda w: config["treatments"]["forester_1"][w.drug]["sigs_names"],
        A_name=lambda w: config["treatments"]["forester_1"][w.drug]["A_name"],
        B_name=lambda w: config["treatments"]["forester_1"][w.drug]["B_name"],
    output:
        ["%s/treatments/forester/{tumor_type}/{drug}_{tumor_type}_%s.pdf" % (R_FOLDER, x) for x in ["clinical", "molecular"]]
    shell:
        """
        Rscript workflow/scripts/02.3_trt_forester.R \
            --cohort {params.cohort} \
            --combination {params.combination} \
            --plot_mode {params.plot_mode} \
            --drug_classes {params.drug_classes} \
            --drug_names {params.drug_names} \
            --drug_table {input.drug_table} \
            --tumor_types {wildcards.tumor_type} \
            --sigs_table {input.sigs_table} \
            --sigs_names {params.sigs_names} \
            --alts_table {input.alts_table} \
            --A_name {params.A_name} \
            --B_name {params.B_name} \
            --outputs {output} \
            --log {log}
        """


rule trt_forester_2:
    log:
        "%s/trt_forester_2_{drug}_{tumor_type}.log" % L_FOLDER
    input:
        drug_table="%s/%s" % (D_FOLDER, FILEPATHS["resources"]["drug_table"]),
        sigs_table=config["data"]["sigs_table"],
        alts_table=config["data"]["alts_table"],
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        combination="union",
        plot_mode=2,
        drug_classes=lambda w: config["treatments"]["forester_2"][w.drug]["drug_classes"],
        drug_names=lambda w: config["treatments"]["forester_2"][w.drug]["drug_names"],
        sigs_names=lambda w: config["treatments"]["forester_2"][w.drug]["sigs_names"],
        A_name=lambda w: config["treatments"]["forester_2"][w.drug]["A_name"],
        B_name=lambda w: config["treatments"]["forester_2"][w.drug]["B_name"],
    output:
        "%s/treatments/forester/{tumor_type}/{drug}_{tumor_type}.pdf" % R_FOLDER
    shell:
        """
        Rscript workflow/scripts/02.3_trt_forester.R \
            --cohort {params.cohort} \
            --combination {params.combination} \
            --plot_mode {params.plot_mode} \
            --drug_classes {params.drug_classes} \
            --drug_names {params.drug_names} \
            --drug_table {input.drug_table} \
            --tumor_types {wildcards.tumor_type} \
            --sigs_table {input.sigs_table} \
            --sigs_names {params.sigs_names} \
            --alts_table {input.alts_table} \
            --A_name {params.A_name} \
            --B_name {params.B_name} \
            --outputs {output} \
            --log {log}
        """


rule trt_table_extra:
    log:
        "%s/trt_table_extra_{col_row}.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort=config["treatments"]["cohort"],
        min_count_cols=config["tumor_types"]["min_count"],
        min_count_rows=25,
        min_frac_rows=0.5
    output:
        "%s/treatments/heatmaps/table_extra_drugs_{col_row}_all_tt.pdf" % R_FOLDER,
        "%s/treatments/heatmaps/table_extra_drugs_{col_row}_sel_tt.pdf" % R_FOLDER,
        "%s/F1d_{col_row}.pdf" % F_FOLDER
    shell:
        """
        Rscript workflow/scripts/02.5_trt_table_extra.R \
            --cohort {params.cohort} \
            --min_count_rows {params.min_count_rows} \
            --min_frac_rows {params.min_frac_rows} \
            --min_count_cols {params.min_count_cols} \
            --col_row {wildcards.col_row} \
            --outputs {output[0]} {output[1]} \
            --output_paper {output[2]} \
            --log {log}
        """
