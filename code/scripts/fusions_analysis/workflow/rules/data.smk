rule aggregate_tables_samples:
    log:
        "workflow/logs/aggregate_tables_samples_{cohort}_{algo}.log"
    input:
        "../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        agg="%s/{cohort}/rna/{algo}/{cohort}_{algo}.tsv.gz" % D_FOLDER
    params:
        output_list="%s/{cohort}/rna/{algo}/sample_list.tsv" % D_FOLDER,
        data_folder=D_FOLDER,
    resources:
        mem_mb=16000,
        partition="shortq",
        time_min=30
    threads: 1
    shell:
        """python workflow/scripts/00.1_aggregate_tables_samples.py \
            --cohort {wildcards.cohort} \
            --algo_folder {wildcards.algo} \
            --data_folder {params.data_folder} \
            --output_list {params.output_list} \
            --output_agg {output.agg} &> {log}
        """


rule aggregate_tables_callers:
    log:
        "workflow/logs/aggregate_tables_callers_{cohort}.log"
    input:
        agg=get_input_aggregate_tables_callers,
        env="../common/logs/setup_conda.done"
    conda: config["setup"]["MetaPrism"]
    output:
        "%s/{cohort}/rna/fusions/{cohort}_aggregated_callers.tsv.gz" % D_FOLDER
    params:
        algos=lambda w: config["data"]["aggregate"][w.cohort]["algos"]
    resources:
        mem_mb=30000,
        partition="shortq",
        time_min=60
    threads: 1
    shell:
        """Rscript workflow/scripts/00.2_aggregate_tables_callers.R \
            --cohort {wildcards.cohort} \
            --algos {params.algos} \
            --output {output} \
            --log {log}"""


rule annotate_fusions_FusionAnnotator_1:
    log:
        "workflow/logs/annotate_fusions_FusionAnnotator_1_{cohort}.log"
    input:
        table="%s/{cohort}/rna/fusions/{cohort}_aggregated_callers.tsv.gz" % D_FOLDER,
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        temp("%s/{cohort}/rna/fusions/{cohort}_aggregated_FusionAnnotator_1.tsv" % D_FOLDER)
    resources:
        mem_mb=16000,
        partition="shortq",
        time_min=60
    threads: 1
    shell:
        """python workflow/scripts/00.3_annotate_fusions_FusionAnnotator_1.py --input {input.table} \
            --output {output} &> {log}"""


rule annotate_fusions_FusionAnnotator_2:
    benchmark:
        "workflow/benchmarks/annotate_fusions_FusionAnnotator_2_{cohort}.tsv"
    log:
        "workflow/logs/annotate_fusions_FusionAnnotator_2_{cohort}.log"
    input:
        "%s/{cohort}/rna/fusions/{cohort}_aggregated_FusionAnnotator_1.tsv" % D_FOLDER
    conda:
        "../envs/FusionAnnotator.yaml"
    output:
        temp("%s/{cohort}/rna/fusions/{cohort}_aggregated_FusionAnnotator_2.tsv" % D_FOLDER)
    params:
        app=config["data"]["fusion_annotator"],
        genome_lib_dir=config["data"]["resources"]["genome_lib_dir"]
    resources:
        mem_mb=16000,
        partition="shortq",
        time_min=60
    threads: 1
    shell:
        """{params.app} --genome_lib_dir {params.genome_lib_dir} \
            --annotate {input} \
            --fusion_name_col Fusion_Id 1> {output} 2> {log}"""


rule annotate_fusions_FusionAnnotator_3:
    log:
        "workflow/logs/annotate_fusions_FusionAnnotator_3_{cohort}.log"
    input:
        fusions="%s/{cohort}/rna/fusions/{cohort}_aggregated_callers.tsv.gz" % D_FOLDER,
        annots="%s/{cohort}/rna/fusions/{cohort}_aggregated_FusionAnnotator_2.tsv" % D_FOLDER,
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        "%s/{cohort}/rna/fusions/{cohort}_annotated_FusionAnnotator.tsv.gz" % D_FOLDER
    resources:
        mem_mb=16000,
        partition="shortq",
        time_min=60
    threads: 1
    shell:
        """python workflow/scripts/00.3_annotate_fusions_FusionAnnotator_3.py --input_fusions {input.fusions} \
            --input_annots {input.annots} \
            --output {output} &> {log}"""


rule annotate_fusions_custom:
    log:
        "workflow/logs/annotate_fusions_custom_{cohort}.log"
    input:
        fusions="%s/{cohort}/rna/fusions/{cohort}_annotated_FusionAnnotator.tsv.gz" % D_FOLDER,
        env="../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        "%s/{cohort}/rna/fusions/{cohort}_annotated.tsv.gz" % D_FOLDER
    params:
        drivers=config["data"]["resources"]["drivers"],
        gencode=config["data"]["resources"]["gencode"],
        fusions_lists=config["data"]["resources"]["fusions_lists"]
    resources:
        mem_mb=8000,
        partition="shortq",
        time_min=15
    threads: 1
    shell:
        """Rscript workflow/scripts/00.4_annotate_fusions_custom.R \
            --input {input.fusions} \
            --fusions_lists {params.fusions_lists} \
            --gencode {params.gencode} \
            --drivers {params.drivers} \
            --output {output} \
            --log {log}"""


rule filter_fusions:
    log:
        "workflow/logs/filter_fusions_{cohort}.log"
    input:
        "%s/{cohort}/rna/fusions/{cohort}_annotated.tsv.gz" % D_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    output:
        fus_filters="%s/{cohort}/rna/fusions/{cohort}_filters.tsv.gz" % D_FOLDER,
        fus="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered.tsv.gz" % D_FOLDER,
        sam="%s/{cohort}/rna/fusions/sample_list.tsv" % D_FOLDER
    resources:
        mem_mb=4000,
        partition="shortq",
        time_min=30
    threads: 1
    shell:
        """python workflow/scripts/00.7_filter_fusions.py \
            --cohort {wildcards.cohort} \
            --input {input} \
            --output_sam {output.sam} \
            --output_fus_filters {output.fus_filters} \
            --output_fus {output.fus} &> {log}
        """


checkpoint oncokb_preprocess:
    input:
        fus="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered.tsv.gz" % D_FOLDER,
        sam="%s/{cohort}/rna/fusions/sample_list.tsv" % D_FOLDER,
        bio="%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER
    output:
        temp(directory("%s/{cohort}/rna/fusions/oncokb_pre" % D_FOLDER))
    benchmark:
        "%s/oncokb_preprocess_{cohort}.tsv" % B_FOLDER
    log:
        "%s/oncokb_preprocess_{cohort}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=10
    shell:
        """
        python workflow/scripts/00.8.1_oncokb_preprocess.py \
            --table_fus {input.fus} \
            --table_sam {input.sam} \
            --table_bio {input.bio} \
            --output {output} &> {log}
        """


rule oncokb_annotate:
    input:
        "%s/{cohort}/rna/fusions/oncokb_pre/{sample}.tsv" % D_FOLDER
    output:
        temp("%s/{cohort}/rna/fusions/oncokb/{sample}.tsv" % D_FOLDER)
    benchmark:
        "%s/oncokb_annotate_{cohort}_{sample}.tsv" % B_FOLDER
    log:
        "%s/oncokb_annotate_{cohort}_{sample}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        code_dir=config["params"]["oncokb"]["code_dir"],
        token=config["params"]["oncokb"]["token"],
        tumor_type=get_tumor_type_mskcc_oncotree
    threads: 1
    resources:
        queue="shortq",
        mem_mb=1000,
        time_min=10
    shell:
        """
        python {params.code_dir}/FusionAnnotator.py \
            -i {input} \
            -b {params.token} \
            -t {params.tumor_type} \
            -o {output} &> {log}
        """


def oncokb_postprocess_input(wildcards):
    checkpoint_output = checkpoints.oncokb_preprocess.get(**wildcards).output[0]
    return expand("%s/{{cohort}}/rna/fusions/oncokb/{sample}.tsv" % D_FOLDER,
                  sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}.tsv")).sample)


rule oncokb_postprocess:
    input:
        fus="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered.tsv.gz" % D_FOLDER,
        sam="%s/{cohort}/rna/fusions/sample_list.tsv" % D_FOLDER,
        bio="%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER,
        okb=oncokb_postprocess_input
    output:
        "%s/{cohort}/rna/fusions/{cohort}_annotated_filtered_oncokb.tsv.gz" % D_FOLDER
    benchmark:
        "%s/oncokb_postprocess_{cohort}.tsv" % B_FOLDER
    log:
        "%s/oncokb_postprocess_{cohort}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=10
    shell:
        """
        python workflow/scripts/00.8.2_oncokb_postprocess.py \
            --table_fus {input.fus} \
            --table_sam {input.sam} \
            --table_bio {input.bio} \
            --oncokb {input.okb} \
            --output {output} &> {log}
        """


checkpoint civic_preprocess:
    input:
        fus="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered.tsv.gz" % D_FOLDER,
        sam="%s/{cohort}/rna/fusions/sample_list.tsv" % D_FOLDER,
        bio="%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER
    output:
        directory("%s/{cohort}/rna/fusions/civic_pre" % D_FOLDER)
    benchmark:
        "%s/civic_preprocess_{cohort}.tsv" % B_FOLDER
    log:
        "%s/civic_preprocess_{cohort}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        python workflow/scripts/00.9.1_civic_preprocess.py \
            --table_fus {input.fus} \
            --table_sam {input.sam} \
            --table_bio {input.bio} \
            --output {output} &> {log}
        """


rule civic_annotate:
    input:
        "%s/{cohort}/rna/fusions/civic_pre/{sample}.tsv" % D_FOLDER
    output:
        "%s/{cohort}/rna/fusions/civic/{sample}.tsv" % D_FOLDER
    benchmark:
        "%s/civic_annotate_{cohort}_{sample}.tsv" % B_FOLDER
    log:
        "%s/civic_annotate_{cohort}_{sample}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    params:
        code_dir=config["params"]["civic"]["code_dir"],
        civic=config["params"]["civic"]["database"],
        rules=config["params"]["civic"]["rules_clean"],
        tumor_type=get_tumor_type_civic
    threads: 1
    resources:
        queue="shortq",
        mem_mb=1000,
        time_min=10
    shell:
        """
        python {params.code_dir}/civic_annotator.py \
            --input {input} \
            --civic {params.civic} \
            --rules {params.rules} \
            --category fus \
            --tumor_types "{params.tumor_type}" \
            --output {output} &> {log}
        """


def civic_postprocess_input(wildcards):
    checkpoint_output = checkpoints.civic_preprocess.get(**wildcards).output[0]
    return expand("%s/{{cohort}}/rna/fusions/civic/{sample}.tsv" % D_FOLDER,
                  sample=glob_wildcards(os.path.join(checkpoint_output, "{sample}.tsv")).sample)


rule civic_postprocess:
    input:
        fus="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered.tsv.gz" % D_FOLDER,
        sam="%s/{cohort}/rna/fusions/sample_list.tsv" % D_FOLDER,
        bio="%s/{cohort}/clinical/curated/bio_{cohort}_in_design_curated.tsv" % D_FOLDER,
        civ=civic_postprocess_input
    output:
        "%s/{cohort}/rna/fusions/{cohort}_annotated_filtered_civic.tsv.gz" % D_FOLDER
    benchmark:
        "%s/civic_postprocess_{cohort}.tsv" % B_FOLDER
    log:
        "%s/civic_postprocess_{cohort}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        python workflow/scripts/00.9.2_civic_postprocess.py \
            --table_fus {input.fus} \
            --table_sam {input.sam} \
            --table_bio {input.bio} \
            --civic {input.civ} \
            --output {output} &> {log}
        """


rule union_ann:
    input:
        civ="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered_civic.tsv.gz" % D_FOLDER,
        okb="%s/{cohort}/rna/fusions/{cohort}_annotated_filtered_oncokb.tsv.gz" % D_FOLDER,
    output:
        "%s/{cohort}/rna/fusions/{cohort}_annotated_filtered_union_ann.tsv.gz" % D_FOLDER
    benchmark:
        "%s/union_ann_{cohort}.tsv" % B_FOLDER
    log:
        "%s/union_ann_{cohort}.log" % L_FOLDER
    conda:
        config["setup"]["MetaPrism"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        python workflow/scripts/00.10_concatenate_fus_annotations.py \
            --civ {input.civ} \
            --okb {input.okb} \
            --output {output} &> {log}
        """
