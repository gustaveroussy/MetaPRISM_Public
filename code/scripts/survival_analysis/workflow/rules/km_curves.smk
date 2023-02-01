rule km_curves_trt_met_rmh_grim:
    wildcard_constraints:
        tumor_type = "|".join([re.escape(x) for x in config["km_curves"]["tumor_types"] if x!="All"]),
    log:
        "%s/km_curves_trt_met_rmh_grim_{tumor_type}.log" % L_FOLDER
    input:
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = config["km_curves"]["cohort"],
        max_strata = 4,
        widths = [10,10,10,10],
        heights = [6,6,6,6]
    output:
        ["%s/km_curves/{tumor_type}/trt_met_rmh_grim/km_curves_%s_{tumor_type}.pdf" % (R_FOLDER, x) for x in ["trt","met","rmh","grim"]]
    shell:
        """
        Rscript workflow/scripts/02.1_km_curves_trt_met_rmh_grim.R \
            --cohort {params.cohort} \
            --tumor_types {wildcards.tumor_type} \
            --max_strata {params.max_strata} \
            --widths {params.widths} \
            --heights {params.heights} \
            --outputs {output} \
            --log {log}
        """


rule km_curves_trt_met_rmh_grim_all:
    log:
        "%s/km_curves_trt_met_rmh_grim_All.log" % L_FOLDER
    input:
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = config["km_curves"]["cohort"],
        max_strata = 4,
        widths = [4.25,4.25,4.25,4.25,5],
        heights = [2.5,2.5,2.5,2.5,4]
    output:
        ["%s/km_curves/All/trt_met_rmh_grim/km_curves_%s_All.pdf" % (R_FOLDER, x) for x in ["trt","met", "rmh","grim"]],
        "%s/F1e.svg" % F_FOLDER
    shell:
        """
        Rscript workflow/scripts/02.1_km_curves_trt_met_rmh_grim.R \
            --cohort {params.cohort} \
            --tumor_types All \
            --max_strata {params.max_strata} \
            --widths {params.widths} \
            --heights {params.heights} \
            --outputs {output} \
            --log {log}
        """

rule km_curves_genome_events:
    wildcard_constraints:
        events = "|".join([re.escape(x) for x in config["km_curves"]["events"]])
    log:
        "%s/km_curves_genome_events_{tumor_type}_{events}.log" % L_FOLDER
    input:
        event_counts = "../../../results/combined_alterations/count/count_by_gene_DNA_RNA_annotated_prism.tsv",
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = config["km_curves"]["cohort"],
        min_counts = lambda w: 20 if w.tumor_type=="dna_cohort" else 12,
        sample_types = ["DNA_T", "RNA_T"],
        width = 5,
        height = 2.5
    output:
        directory("%s/km_curves/{tumor_type}/by_{events}" % R_FOLDER)
    shell:
        """
        mkdir -p {output}
        Rscript workflow/scripts/02.2_km_curves_genome_events.R \
            --cohort {params.cohort} \
            --tumor_types {wildcards.tumor_type} \
            --event_counts {input.event_counts} \
            --sample_types {params.sample_types} \
            --min_counts {params.min_counts} \
            --granularity {wildcards.events} \
            --width {params.width} \
            --height {params.height} \
            --output {output} \
            --log {log}
        """


rule km_curves_tme_subtypes:
    log:
        "%s/km_curves_tme_subtypes_{tumor_type}.log" % L_FOLDER
    input:
        tme_subtypes = config["data"]["prism"]["rna_gex_tme_bagaev"],
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    params:
        cohort = config["km_curves"]["cohort"],
        width = 3.25,
        height = 2.25
    output:
        "%s/km_curves/{tumor_type}/by_subtype/km_curves_by_tme_subtypes.pdf" % R_FOLDER
    shell:
        """
        Rscript workflow/scripts/02.3_km_curves_tme_subtypes.R \
            --cohort {params.cohort} \
            --tumor_types {wildcards.tumor_type} \
            --tme_subtypes {input.tme_subtypes} \
            --width {params.width} \
            --height {params.height} \
            --output {output} \
            --log {log}
        """
