rule sankey_plot_levels:
    log:
        "%s/sankey_plot_levels_{cohort}.log" % L_FOLDER
    input:
        cln_table = lambda w: config["data"]["cln"][w.cohort],
        alt_table = "%s/alterations/aggregated_alterations_{cohort}.tsv" % R_FOLDER,
        sam_table = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        "%s/other_plots/sankey_levels_{cohort}.svg" % R_FOLDER,
    params:
        min_count = 10,
        width = 950,
        height = 1250,
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        python workflow/scripts/03.1_sankey_plot_levels.py \
            --cln_table {input.cln_table} \
            --alt_table {input.alt_table} \
            --sam_table {input.sam_table} \
            --min_count {params.min_count} \
            --width {params.width} \
            --height {params.height} \
            --output {output} &> {log}
        """


rule violins_drivers_dna:
    log:
        "%s/violins_drivers_dna.log" % L_FOLDER
    input:
        counts_tt = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        alt_tables = expand("%s/alterations/aggregated_alterations_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"])
    conda:
        config["setup"]["MetaPrism"]
    output:
        out = "%s/other_plots/violins_drivers_dna.eps" % R_FOLDER,
        out_paper = "%s/F3b.eps" % F_FOLDER
    params:
        cohorts = config["data"]["cohorts"]
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        Rscript workflow/scripts/03.2_violins_drivers_dna.R \
            --cohorts {params.cohorts} \
            --counts_tt {input.counts_tt} \
            --alt_tables {input.alt_tables} \
            --output {output.out} \
            --output_paper {output.out_paper} \
            --log {log}
        """


rule barplot_multihits:
    log:
        "%s/barplot_multihits_{gene_type}.log" % L_FOLDER
    input:
        alt_table = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        cna_tables = expand("%s/{cohort}/wes/somatic_cna/somatic_calls.tsv.gz" % D_FOLDER,
            cohort=config["data"]["cohorts"]),
        gen_table = "%s/resources/curated/cancer_genes_curated.tsv" % D_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    output:
        table="%s/other_plots/table_multihits_{gene_type}.tsv" % R_FOLDER,
        plots=["%s/other_plots/barplot_multihits_{gene_type}_%s.pdf" % (R_FOLDER,x) \
            for x in ["All", "BLCA", "BRCA", "LUAD", "PAAD", "PRAD"]]
        papers=["%s/FS12_13_{gene_type}_%s.pdf" % (F_FOLDER, x) \
            for x in ["All", "BLCA", "BRCA", "LUAD", "PAAD", "PRAD"]]
    params:
        cohorts = config["data"]["cohorts"],
        output_width=650,
        output_height=250
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    shell:
        """
        Rscript workflow/scripts/03.3_barplot_multihits.R \
            --cohorts {params.cohorts} \
            --alt_table {input.alt_table} \
            --cna_tables {input.cna_tables} \
            --gen_table {input.gen_table} \
            --gene_type {wildcards.gene_type} \
            --output_width {params.output_width} \
            --output_height {params.output_height} \
            --output_table {output.table} \
            --output_plots {output.plots} \
            --output_papers {output.papers} \
            --log {log}
        """
