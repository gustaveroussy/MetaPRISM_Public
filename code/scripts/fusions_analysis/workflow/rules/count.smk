rule count_total:
    log:
        "workflow/logs/count_total_for_all_{cohort}.log"
    input:
        samples = lambda w: "%s/%s/rna/fusions/sample_list.tsv" % (D_FOLDER, w.cohort),
        fusions = lambda w: "%s/%s" % (D_FOLDER, FILEPATHS[w.cohort]["rna_fus"]["fusions"]),
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/count/count_total_RNA_all_{cohort}.tsv" % R_FOLDER
    threads: 1
    shell:
        """
        Rscript workflow/scripts/02.1_count_total.R \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --fusions {input.fusions} \
            --output {output} \
            --log {log}
        """
