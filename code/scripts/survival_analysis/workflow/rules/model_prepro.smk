rule gather_all_features:
    log:
        "%s/gather_all_features_{cohort}.log" % L_FOLDER
    input:
        cln = lambda w: config["data"][w.cohort]["cln"],
        cln_brca = lambda w: config["data"][w.cohort]["cln_brca"],
        cln_covariates = lambda w: config["data"][w.cohort]["cln_covariates"],
        driver_genes = lambda w: config["data"][w.cohort]["driver_genes"],
        dna_cna_chr_arm = lambda w: config["data"][w.cohort]["dna_cna_chr_arm"],
        dna_cna_per_gene = lambda w: config["data"][w.cohort]["dna_cna_per_gene"],
        dna_cna_summary_stats = lambda w: config["data"][w.cohort]["dna_cna_summary_stats"],
        dna_mut_signatures = lambda w: config["data"][w.cohort]["dna_mut_signatures"],
        dna_mut_counts_gene = lambda w: config["data"][w.cohort]["dna_mut_counts_gene"],
        dna_mut_counts_pathway = lambda w: config["data"][w.cohort]["dna_mut_counts_pathway"],
        dna_mut_counts_total = lambda w: config["data"][w.cohort]["dna_mut_counts_total"],
        dna_alt_counts_gene = lambda w: config["data"][w.cohort]["dna_alt_counts_gene"],
        dna_alt_counts_pathway = lambda w: config["data"][w.cohort]["dna_alt_counts_pathway"],
        dna_alt_counts_total = lambda w: config["data"][w.cohort]["dna_alt_counts_total"],
        dna_rna_alt_counts_gene = lambda w: config["data"][w.cohort]["dna_rna_alt_counts_gene"],
        dna_rna_alt_counts_pathway = lambda w: config["data"][w.cohort]["dna_rna_alt_counts_pathway"],
        dna_rna_alt_counts_total = lambda w: config["data"][w.cohort]["dna_rna_alt_counts_total"],
        dna_msi_status_score = lambda w: config["data"][w.cohort]["dna_msi_status_score"],
        rna_gex_tf_signatures = lambda w: config["data"][w.cohort]["rna_gex_tf_signatures"],
        rna_gex_tme_bagaev = lambda w: config["data"][w.cohort]["rna_gex_tme_bagaev"],
        rna_gex_ssgsea_pathway = lambda w: config["data"][w.cohort]["rna_gex_ssgsea_pathway"],
        rna_fus_counts_total = lambda w: config["data"][w.cohort]["rna_fus_counts_total"],
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        dat = "%s/data_{cohort}/all_features/data.tsv.gz" % R_FOLDER,
        cov = "%s/data_{cohort}/all_features/covs.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "00:45:00"
    shell:
        """
        Rscript ../common/scripts/model_prepro_gather_all_features.R \
            --cohort {wildcards.cohort} \
            --cln {input.cln} \
            --cln_brca {input.cln_brca} \
            --cln_covariates {input.cln_covariates} \
            --driver_genes {input.driver_genes} \
            --dna_cna_chr_arm {input.dna_cna_chr_arm} \
            --dna_cna_per_gene {input.dna_cna_per_gene} \
            --dna_cna_summary_stats {input.dna_cna_summary_stats} \
            --dna_mut_signatures {input.dna_mut_signatures} \
            --dna_mut_counts_gene {input.dna_mut_counts_gene} \
            --dna_mut_counts_pathway {input.dna_mut_counts_pathway} \
            --dna_mut_counts_total {input.dna_mut_counts_total} \
            --dna_alt_counts_gene {input.dna_alt_counts_gene} \
            --dna_alt_counts_pathway {input.dna_alt_counts_pathway} \
            --dna_alt_counts_total {input.dna_alt_counts_total} \
            --dna_rna_alt_counts_gene {input.dna_rna_alt_counts_gene} \
            --dna_rna_alt_counts_pathway {input.dna_rna_alt_counts_pathway} \
            --dna_rna_alt_counts_total {input.dna_rna_alt_counts_total} \
            --dna_msi_status_score {input.dna_msi_status_score} \
            --rna_gex_tf_signatures {input.rna_gex_tf_signatures} \
            --rna_gex_tme_bagaev {input.rna_gex_tme_bagaev} \
            --rna_gex_ssgsea_pathway {input.rna_gex_ssgsea_pathway} \
            --rna_fus_counts_total {input.rna_fus_counts_total} \
            --output_dat {output.dat} \
            --output_cov {output.cov} \
            --log {log}
        """


rule select_sub_features:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES])
    log:
        "%s/select_sub_features_{cohort}_{samples}_{features}.log" % L_FOLDER
    input:
        dat = "%s/data_{cohort}/all_features/data.tsv.gz" % R_FOLDER,
        cov = "%s/data_{cohort}/all_features/covs.tsv.gz" % R_FOLDER,
        sam = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        cnt = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        config_yaml = "config/config.yaml",
    output:
        dat = "%s/data_{cohort}/sub_features/{samples}/{features}/original/data.original.tsv.gz" % R_FOLDER,
        cov = "%s/data_{cohort}/sub_features/{samples}/{features}/original/covs.original.tsv.gz" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "00:15:00"
    shell:
        """
        Rscript ../common/scripts/model_prepro_select_sub_features.R  \
            --input_dat {input.dat} \
            --input_cov {input.cov} \
            --input_sam {input.sam} \
            --input_cnt {input.cnt} \
            --config_yaml {params.config_yaml} \
            --config_section models \
            --features {wildcards.features} \
            --samples {wildcards.samples} \
            --output_dat {output.dat} \
            --output_cov {output.cov} \
            --log {log}
        """


checkpoint group_encode_check_impute:
    wildcard_constraints:
        cohort = "|".join([re.escape(x) for x in COHORTS]),
        samples = "|".join([re.escape(x) for x in SAMPLES]),
        features = "|".join([re.escape(x) for x in FEATURES])
    benchmark:
        "%s/group_encode_check_impute_{cohort}_{samples}_{features}.log" % B_FOLDER
    log:
        "%s/group_encode_check_impute_{cohort}_{samples}_{features}.log" % L_FOLDER
    input:
        dat = "%s/data_{cohort}/sub_features/{samples}/{features}/original/data.original.tsv.gz" % R_FOLDER,
        cov = "%s/data_{cohort}/sub_features/{samples}/{features}/original/covs.original.tsv.gz" % R_FOLDER,
    conda:
        config["setup"]["MetaPrism"]
    params:
        spread_discrete_ordered = config["models"]["preprocess"]["spread_discrete_ordered"],
        min_size = config["models"]["preprocess"]["min_size"],
        seed = config["models"]["preprocess"]["seed"],
        n_imputations = config["models"]["preprocess"]["n_imputations"],
        n_iterations = config["models"]["preprocess"]["n_iterations"],
    threads: 4
    output:
        dir = directory("%s/data_{cohort}/sub_features/{samples}/{features}/processed" % R_FOLDER),
        cov_final = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.final.tsv.gz" % R_FOLDER,
        cov_removed = "%s/data_{cohort}/sub_features/{samples}/{features}/processed/covs.removed.tsv.gz" % R_FOLDER,
        pmc = "%s/data_{cohort}/sub_features/{samples}/{features}/original/plot_missing_data_counts.pdf" % R_FOLDER,
        pmp = "%s/data_{cohort}/sub_features/{samples}/{features}/original/plot_missing_data_patterns.pdf" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "00:30:00"
    shell:
        """
        Rscript ../common/scripts/model_prepro_group_encode_check_impute.R \
            --input_dat {input.dat} \
            --input_cov {input.cov} \
            --spread_discrete_ordered {params.spread_discrete_ordered} \
            --config_yaml config/config.yaml \
            --config_section models \
            --features {wildcards.features} \
            --min_size {params.min_size} \
            --seed {params.seed} \
            --n_imputations {params.n_imputations} \
            --n_iterations {params.n_iterations} \
            --n_cores {threads} \
            --output_pmc {output.pmc} \
            --output_pmp {output.pmp} \
            --output_dir {output.dir} \
            --log {log}
        """
