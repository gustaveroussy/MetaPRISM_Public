rule update_gene_symbols:
    log:
        "%s/update_gene_symbols_{gene_table}.log" % L_FOLDER
    input:
        gene = lambda w: config["resources"]["gene_tables"][w.gene_table],
        hgnc = config["resources"]["hgnc"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/gene_tables/{gene_table}_updated.tsv" % R_FOLDER
    params:
        id_input = "Ensembl_Gene_Id",
        id_biomart = "ensembl_gene_id"
    resources:
        partition = "cpu_short",
        mem_mb = 2000,
        time = "00:30:00"
    threads: 1
    shell:
        """
        python workflow/scripts/01.1_update_gene_symbols.py \
            --gene {input.gene} \
            --hgnc {input.hgnc} \
            --gene_name Hugo_Symbol \
            --gene_ensembl_id Ensembl_Gene_Id \
            --output {output} &> {log}
        """


rule signatures_scores_mfp_data:
    benchmark:
        "%s/signatures_scores_mfp_{cohort}.tsv" % B_FOLDER
    log:
        "%s/signatures_scores_mfp_{cohort}.log" % L_FOLDER
    input:
        expression_data = get_expression_data,
        expression_summary = get_expression_summary,
        bio_data = "%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER,
        gene_table = get_gene_table_udpated,
        gene_selection = get_gene_selection_udpated,
        signatures = config["signatures"]["gene_signatures"]["bagaev_2021"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/{cohort}/tables/signatures_bagaev_2021_mfp.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = lambda w: config["signatures"]["mem_mb"][w.cohort],
        time = "01:00:00"
    threads: 1
    shell:
        """
        python workflow/scripts/02.1_signatures_scores_mfp.py \
            --expression_data {input.expression_data} \
            --expression_summary {input.expression_summary} \
            --bio_data {input.bio_data} \
            --gene_table {input.gene_table} \
            --gene_selection {input.gene_selection} \
            --signatures {input.signatures} \
            --output {output} &> {log}
        """


rule signatures_scores_gsva_data:
    benchmark:
        "%s/signatures_scores_gsva_{cohort}.tsv" % B_FOLDER
    input:
        expression_data = get_expression_data,
        expression_summary = get_expression_summary,
        bio_data = "%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER,
        gene_table = get_gene_table_udpated,
        gene_selection = get_gene_selection_udpated,
        signatures = config["signatures"]["gene_signatures"]["bagaev_2021"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/signatures_scores_gsva_{cohort}.log" % L_FOLDER
    output:
        "%s/{cohort}/tables/signatures_bagaev_2021_gsva.tsv" % R_FOLDER
    params:
        gsva_method = config["signatures"]["gsva_method"],
        gsva_ssgsea_norm = config["signatures"]["gsva_ssgsea_norm"],
        gsva_tau = config["signatures"]["gsva_tau"]
    resources:
        partition = "cpu_med",
        mem_mb = lambda w: config["signatures"]["mem_mb"][w.cohort],
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.1_signatures_scores_gsva.R \
            --expression_data {input.expression_data} \
            --expression_summary {input.expression_summary} \
            --bio_data {input.bio_data} \
            --gene_table {input.gene_table} \
            --gene_selection {input.gene_selection} \
            --signatures {input.signatures} \
            --gsva_method {params.gsva_method} \
            --gsva_ssgsea_norm {params.gsva_ssgsea_norm} \
            --gsva_tau {params.gsva_tau} \
            --output {output} \
            --log {log}
        """
