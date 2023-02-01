rule double_htmp_treatments:
    log:
        "%s/double_htmp_treatments_{tumor_type}_{level}.log" % L_FOLDER
    input:
        drug_table = "%s/resources/drug_tables/Table_Drugs_v7.xlsx" % D_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = "prism"
    threads: 1
    output:
        "%s/plots/double_htmp_{tumor_type}_{level}.pdf" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 1000,
        time = "00:15:00"
    shell:
        """
        python workflow/scripts/01.1_double_htmp_treatments.py \
            --cohort {params.cohort} \
            --drug_table {input.drug_table} \
            --tumor_type {wildcards.tumor_type} \
            --level {wildcards.level} \
            --output {output} &> {log}
        """


rule process_annotations:
    log:
        "%s/process_annotations.log" % L_FOLDER
    input:
        alterations_before_annot = "%s/annotations/aggregated_alterations_prism_072022.xlsx"  % R_FOLDER,
        alterations_after_annot = "%s/annotations/aggregated_alterations_prism_annotated_LF_072022.xlsx"  % R_FOLDER,
        alterations_current = "../../../results/combined_alterations/alterations/aggregated_alterations_prism_all.tsv",
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    output:
        alterations = "%s/annotations/aggregated_alterations_prism_all_annot.tsv" % R_FOLDER,
        rules_annot = "%s/annotations/rules_annotation_LF.xlsx" % R_FOLDER,
        treatments = "%s/annotations/treatment_resistances_prism.xlsx" % R_FOLDER,
    resources:
        partition = "cpu_med",
        mem_mb = 1000,
        time = "00:15:00"
    shell:
        """
        python workflow/scripts/01.2_process_annotations.py \
            --alterations_before_annot {input.alterations_before_annot} \
            --alterations_after_annot {input.alterations_after_annot} \
            --alterations_current {input.alterations_current} \
            --output_alterations {output.alterations} \
            --output_rules_annot {output.rules_annot} \
            --output_treatments {output.treatments} &>  {log}
        """


rule oncoplot_treatments:
    log:
        "%s/oncoplot_treatments_{tumor_type}_{level}.log" % L_FOLDER
    input:
        alterations = "%s/annotations/aggregated_alterations_prism_all_annot.tsv" % R_FOLDER,
        drug_table = "%s/resources/drug_tables/Table_Drugs_v7.xlsx" % D_FOLDER,
        drug_rules = "%s/resources/drug_tables/Table_Drugs_Groups_Oncoplot_29072022.xlsx" % D_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = "prism",
        threshold_drug = 5,
        threshold_odds = 4,
        threshold_stks = 5
    threads: 1
    output:
        table = "%s/plots/comuts/table_{tumor_type}_{level}.xlsx" % R_FOLDER,
        plot = "%s/plots/comuts/oncoplot_{tumor_type}_{level}.pdf" % R_FOLDER,
        plot_paper = "%s/F6_{tumor_type}_{level}.eps" % F_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 1000,
        time = "00:15:00"
    shell:
        """
        python workflow/scripts/01.3_oncoplot_treatments.py \
            --cohort {params.cohort} \
            --alterations {input.alterations} \
            --drug_table {input.drug_table} \
            --drug_rules {input.drug_rules} \
            --tumor_type {wildcards.tumor_type} \
            --level {wildcards.level} \
            --threshold_drug {params.threshold_drug} \
            --threshold_odds {params.threshold_odds} \
            --threshold_stks {params.threshold_stks} \
            --output_table {output.table} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} &> {log}
        """
