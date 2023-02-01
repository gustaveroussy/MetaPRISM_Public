rule count_by_gene:
    log:
        "workflow/logs/count_by_gene_{data_type}_annotated_{cohort}.log"
    input:
        alterations = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_gene_{data_type}_annotated_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/04.1_count_by_gene.R \
            --alterations {input.alterations} \
            --samples {input.samples} \
            --data_type {wildcards.data_type} \
            --cohort {wildcards.cohort} \
            --output {output} \
            --log {log}
        """


rule count_by_pathway:
    wildcard_constraints:
        pathway = "|".join([re.escape(x) for x in config["count"]["pathways"]])
    log:
        "workflow/logs/count_by_pathway_{pathway}_{data_type}_annotated_{cohort}.log"
    input:
        alterations = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        pathways = lambda w: config["count"]["pathways"][w.pathway],
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_pathway_{pathway}_{data_type}_annotated_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/04.2_count_by_pathway.R \
            --alterations {input.alterations} \
            --samples {input.samples} \
            --data_type {wildcards.data_type} \
            --cohort {wildcards.cohort} \
            --pathways {input.pathways} \
            --output {output} \
            --log {log}
        """


rule count_total:
    log:
        "workflow/logs/count_total_{data_type}_annotated_{cohort}.log"
    input:
        alterations = "%s/alterations/aggregated_alterations.tsv" % R_FOLDER,
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_total_{data_type}_annotated_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 8000,
        time = "01:00:00"
    threads: 1
    shell:
        """
        Rscript workflow/scripts/04.3_count_total.R \
            --alterations {input.alterations} \
            --samples {input.samples} \
            --data_type {wildcards.data_type} \
            --cohort {wildcards.cohort} \
            --output {output} \
            --log {log}
        """
