rule count_by_pathway:
    wildcard_constraints:
        cna_mode = "|".join([re.escape(x) for x in config["count"]["cna_modes"]]),
        pathway = "|".join([re.escape(x) for x in config["count"]["pathways"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_by_pathway_{cna_mode}_{pathway}_{cohort}.log"
    input:
        cnas = lambda w: config["data"]["cnas"][w.cna_mode][w.cohort],
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        pathways = lambda w: config["data"]["pathways"][w.pathway],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_pathway_{cna_mode}_{pathway}_{cohort}.tsv" % R_FOLDER
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.1_count_by_pathway.R \
            --cnas {input.cnas} \
            --samples {input.samples} \
            --pathways {input.pathways} \
            --output {output} \
            --log {log}
        """


rule count_by_gene:
    wildcard_constraints:
        cna_mode = "|".join([re.escape(x) for x in config["count"]["cna_modes"]]),
        pathway = "|".join([re.escape(x) for x in config["count"]["genes"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_by_pathway_{cna_mode}_{pathway}_{cohort}.log"
    input:
        cnas = lambda w: config["data"]["cnas"][w.cna_mode][w.cohort],
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        genes = lambda w: config["data"]["genes"][w.pathway],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_gene_{cna_mode}_{pathway}_{cohort}.tsv" % R_FOLDER
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.2_count_by_gene.R \
            --cnas {input.cnas} \
            --samples {input.samples} \
            --genes {input.genes} \
            --output {output} \
            --log {log}
        """
