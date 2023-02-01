rule count_by_gene:
    wildcard_constraints:
        selection = "|".join([re.escape(x) for x in config["count"]["selections"]]),
        gene = "|".join([re.escape(x) for x in config["count"]["genes"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_by_gene_{gene}_for_{selection}_{cohort}.log"
    input:
        mutations = get_mutations_file,
        samples = lambda w: config["data"]["mut"]["sam"][w.cohort],
        genes = lambda w: config["data"]["genes"][w.gene],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_gene_{gene}_DNA_{selection}_{cohort}.tsv" % R_FOLDER
    threads: 1
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    shell:
        """
        Rscript workflow/scripts/04.1_count_by_gene.R \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --genes {input.genes} \
            --selection {wildcards.selection} \
            --output {output} \
            --log {log}
        """

rule count_by_pathway:
    wildcard_constraints:
        selection = "|".join([re.escape(x) for x in config["count"]["selections"]]),
        pathway = "|".join([re.escape(x) for x in config["count"]["pathways"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_by_pathway_{pathway}_for_{selection}_{cohort}.log"
    input:
        mutations = get_mutations_file,
        samples = lambda w: config["data"]["mut"]["sam"][w.cohort],
        pathways = lambda w: config["data"]["pathways"][w.pathway],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_by_pathway_{pathway}_DNA_{selection}_{cohort}.tsv" % R_FOLDER
    threads: 1
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    shell:
        """
        Rscript workflow/scripts/04.2_count_by_pathway.R \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --pathways {input.pathways} \
            --selection {wildcards.selection} \
            --output {output} \
            --log {log}
        """


rule count_total:
    wildcard_constraints:
        selection = "|".join([re.escape(x) for x in config["count"]["selections"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_total_for_{selection}_{cohort}.log"
    input:
        mutations = get_mutations_file,
        samples = lambda w: config["data"]["mut"]["sam"][w.cohort],
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_total_DNA_{selection}_{cohort}.tsv" % R_FOLDER
    threads: 1
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    shell:
        """
        Rscript workflow/scripts/04.3_count_total.R \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --selection {wildcards.selection} \
            --output {output} \
            --log {log}
        """
