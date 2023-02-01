rule mutpanning_compiling:
    log:
        "workflow/logs/mutpanning_compiling.log"
    output:
        touch("workflow/logs/mutpanning_compiling.done")
    params:
        app = config["mutpanning"]["app"]
    threads: 1
    shell:
        """
        module load intel/19.0.3/gcc-4.8.5
        module load openjdk/11.0.2/intel-19.0.3.199
        javac -classpath {params.app}/commons-math3-3.6.1.jar:{params.app}/jdistlib-0.4.5-bin.jar {params.app}/*.java
        """


rule mutpanning_inputs:
    log:
        "%s/mutpanning_inputs_{cohort}_{tumor_type}.log" % L_FOLDER
    input:
        mutations = lambda w: config["data"]["mut"]["all"][w.cohort],
        samples = "%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        ann = "%s/mutpanning/{cohort}/{tumor_type}_inputs/samples_{tumor_type}.txt" % R_FOLDER,
        maf = "%s/mutpanning/{cohort}/{tumor_type}_inputs/mutations_{tumor_type}.maf" % R_FOLDER
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.1_mutpanning_inputs.R \
            --cohort {wildcards.cohort} \
            --samples {input.samples} \
            --mutations {input.mutations} \
            --tumor_type {wildcards.tumor_type} \
            --output_ann {output.ann} \
            --output_maf {output.maf} \
            --log {log}
        """


rule mutpanning_run:
    benchmark:
        "workflow/benchmarks/mutpanning_{cohort}_{tumor_type}.tsv"
    log:
        "workflow/logs/mutpanning_{cohort}_{tumor_type}.log"
    input:
        compiling = rules.mutpanning_compiling.output,
        ann = "%s/mutpanning/{cohort}/{tumor_type}_inputs/samples_{tumor_type}.txt" % R_FOLDER,
        maf = "%s/mutpanning/{cohort}/{tumor_type}_inputs/mutations_{tumor_type}.maf" % R_FOLDER
    output:
        touch("%s/mutpanning_{cohort}_{tumor_type}.done" % L_FOLDER)
    params:
        dir = "%s/mutpanning/{cohort}/{tumor_type}" % R_FOLDER,
        app = config["mutpanning"]["app"],
        hg19 = config["mutpanning"]["hg19"]
    resources:
        partition = "cpu_med",
        mem_mb = 16000,
        time = "03:00:00"
    threads: 1
    shell:
        """
        mkdir -p {params.dir}
        module load openjdk/11.0.2/intel-19.0.3.199
        java -Xmx8G \
         -classpath {params.app}/commons-math3-3.6.1.jar:{params.app}/jdistlib-0.4.5-bin.jar:{params.app} MutPanning \
         "{params.dir}/" "{input.maf}" "{input.ann}" "{params.hg19}/"
        """


rule mutpanning_results:
    log:
        "%s/mutpanning_results.log" % L_FOLDER
    input:
        done=expand("%s/mutpanning_{cohort}_{tumor_type}.done" % L_FOLDER, cohort=config["data"]["cohorts"],
            tumor_type=config["mutpanning"]["tumor_types"]),
        counts="%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        sams=expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER, cohort=config["data"]["cohorts"]),
        muts=[config["data"]["mut"]["all"][cohort] for cohort in config["data"]["cohorts"]],
        env="../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        plot="%s/mutpanning/table_extra_mutpanning.pdf" % R_FOLDER,
        plot_paper="%s/FS12.pdf" % F_FOLDER,
        data="%s/mutpanning/results_mutpanning.tsv" % R_FOLDER
    params:
        cohorts=config["data"]["cohorts"],
        mpgs=["%s/mutpanning/%s" % (R_FOLDER, cohort) for cohort in config["data"]["cohorts"]],
    threads: 1
    shell:
        """
        Rscript workflow/scripts/01.2_mutpanning_results.R \
            --cohorts {params.cohorts} \
            --counts {input.counts} \
            --sams {input.sams} \
            --muts {input.muts} \
            --mpgs {params.mpgs} \
            --output_plot {output.plot} \
            --output_plot_paper {output.plot_paper} \
            --output_data {output.data} \
            --log {log}
        """
