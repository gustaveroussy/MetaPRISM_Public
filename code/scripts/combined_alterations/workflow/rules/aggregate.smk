rule aggregate_alterations_across_modalities:
    log:
        "%s/aggregate_alterations_across_modalities_{cohort}.log" % L_FOLDER
    input:
        bio = lambda w: config["data"]["bio"][w.cohort],
        cln = lambda w: config["data"]["cln"][w.cohort],
        cna = lambda w: config["data"]["cna"][w.cohort],
        fus = lambda w: config["data"]["fus"][w.cohort],
        msi = lambda w: config["data"]["msi"][w.cohort],
        mut = lambda w: config["data"]["mut"][w.cohort],
        tmb = lambda w: config["data"]["tmb"][w.cohort],
        exp_arv7 = lambda w: config["data"]["exp_arv7"][w.cohort],
        gen = config["data"]["resources"]["cancer_genes"],
        target_bed = config["data"]["resources"]["target_bed"],
        drug = config["data"]["resources"]["drug"],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        best="%s/alterations/aggregated_alterations_{cohort}.tsv" % R_FOLDER,
        all="%s/alterations/aggregated_alterations_{cohort}_all.tsv" % R_FOLDER,
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.1_aggregate_alterations_across_modalities.R \
            --bio {input.bio} \
            --cln {input.cln} \
            --cna {input.cna} \
            --fus {input.fus} \
            --msi {input.msi} \
            --mut {input.mut} \
            --tmb {input.tmb} \
            --exp_arv7 {input.exp_arv7} \
            --target_bed {input.target_bed} \
            --drug {input.drug} \
            --output_best {output.best} \
            --output_all {output.all} \
            --log {log}
        """


rule aggregate_alterations_across_cohorts:
    log:
        "%s/aggregate_alterations_across_cohorts.log" % L_FOLDER
    input:
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        tables = expand("%s/alterations/aggregated_alterations_{cohort}.tsv" % R_FOLDER,
            cohort=config["data"]["cohorts"]),
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/alterations/aggregated_alterations.tsv" % R_FOLDER
    params:
        cohorts = config["data"]["cohorts"]
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.2_aggregate_alterations_across_cohorts.R \
            --cohorts {params.cohorts} \
            --counts {input.counts} \
            --tables {input.tables} \
            --output {output} \
            --log {log}
        """
