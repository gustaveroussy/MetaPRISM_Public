rule count_total:
    wildcard_constraints:
        selection_mut = "|".join([re.escape(x) for x in config["count"]["selection_muts"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]])
    log:
        "workflow/logs/count_total_for_{selection_mut}_{cohort}.log"
    input:
        mutations = get_mutations_file,
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_total_for_{selection_mut}_{cohort}.tsv" % R_FOLDER
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.1_count_total.R \
            --cohort {wildcards.cohort} \
            --mutations {input.mutations} \
            --samples {input.samples} \
            --selection_mut {wildcards.selection_mut} \
            --output {output} \
            --log {log}
        """
