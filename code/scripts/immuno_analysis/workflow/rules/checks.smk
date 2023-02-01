rule mfp_tcga_check_expression:
    input:
        expression_prism = get_expression_data(pd.Series({"cohort": "tcga"})),
        summary_prism = get_expression_summary(pd.Series({"cohort": "tcga"})),
        gene_table_prism = get_gene_table_udpated(pd.Series({"cohort": "tcga"})),
        expression_xena = config["resources"]["ucsc_xena"]["expression"],
        gene_table_xena = get_gene_table_udpated(pd.Series({"cohort": "ucsc_xena"})),
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_tcga_check_expression.log" % L_FOLDER
    output:
        "%s/tcga/checks/mfp_tcga_check_expression.png" % R_FOLDER
    threads: 1
    resources:
        mem_mb = 48000
    shell:
        """
        python workflow/scripts/03.1_mfp_tcga_check_expression.py \
            --expression_prism {input.expression_prism} \
            --summary_prism {input.summary_prism} \
            --gene_table_prism {input.gene_table_prism} \
            --expression_xena {input.expression_xena} \
            --gene_table_xena {input.gene_table_xena} \
            --output {output} &> {log}
        """


rule mfp_tcga_check_signatures:
    input:
        bagaev_scores = config["resources"]["bagaev_2021"]["signatures"],
        bagaev_annots = config["resources"]["bagaev_2021"]["annotation"],
        prism_scores = "%s/tcga/tables/signatures_mfp_model.tsv" % R_FOLDER,
        prism_annots = "%s/tcga/clinical/curated_other/bio_tcga_all_curated.tsv" % D_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    log:
        "%s/mfp_tcga_check_signatures.log" % L_FOLDER
    output:
        "%s/tcga/checks/mfp_tcga_check_signatures.png" % R_FOLDER
    params:
        method = "pearson"
    threads: 1
    shell:
        """
        python workflow/scripts/03.3_mfp_tcga_check_signatures.py \
            --bagaev_scores {input.bagaev_scores} \
            --bagaev_annots {input.bagaev_annots} \
            --prism_scores {input.prism_scores} \
            --prism_annots {input.prism_annots} \
            --method {params.method} \
            --output {output} &> {log}
        """

