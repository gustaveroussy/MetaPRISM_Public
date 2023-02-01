rule bs_donut_plots:
    log:
        "%s/bs_donut_plots.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohorts=config["tumor_types"]["cohorts"],
        width=405,
        height=405
    output:
        expand("%s/biopsy_sites/donut_plot_tumor_types_{cohort}.svg" % R_FOLDER,cohort=config["tumor_types"]["cohorts"]),
        "%s/F1a_right.svg" % F_FOLDER
    shell:
        """
        python workflow/scripts/01.3_bs_donut_plots.py \
            --cohorts {params.cohorts} \
            --width {params.width} \
            --height {params.height} \
            --outputs {output} &> {log}
        """

rule met_table_extra:
    log:
        "%s/met_table_extra.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort="prism",
        min_count=config["tumor_types"]["min_count"]
    output:
        out=["%s/metastatic_sites/table_extra_met_sites_%s_tt.pdf" % (R_FOLDER, x) for x in ["all", "sel"]],
        paper="%s/FS2.pdf" % F_FOLDER
    shell:
        """
        Rscript workflow/scripts/03.1_met_table_extra.R \
            --cohort {params.cohort} \
            --min_count {params.min_count} \
            --outputs {output} \
            --output_paper {output.paper} \
            --log {log}
        """


rule violins_age_drug_surv_met:
    log:
        "%s/violins_age_drug_surv_met.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort="prism",
        min_count=config["tumor_types"]["min_count"]
    output:
        out=["%s/violins/violins_age_drug_surv_met.eps" % R_FOLDER],
        out_paper=["%s/F1c.eps" % F_FOLDER],
        out_table="%s/violins/violins_medians.tsv" % R_FOLDER
    shell:
        """
        Rscript workflow/scripts/03.2_violins_age_drug_surv_met.R \
            --cohort {params.cohort} \
            --min_count {params.min_count} \
            --outputs {output.out} \
            --outputs_paper {output.out_paper} \
            --output_table {output.out_table} \
            --log {log}
        """


rule venn_molecular:
    log:
        "%s/venn_molecular.log" % L_FOLDER
    input:
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort="prism",
    output:
        out="%s/misc/venn_molecular.pdf" % R_FOLDER
    shell:
        """
        python -u workflow/scripts/03.4_venn_molecular.py \
            --cohort {params.cohort} \
            --output {output.out} \
            --output_paper {output.out_paper} &> {log}
        """
