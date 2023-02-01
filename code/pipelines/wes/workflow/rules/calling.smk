####
#### Germline mutations ####
####

# Call germline SNPs and indels via local re-assembly of haplotypes
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
rule germline_haplotype_caller:
    input:
        ref=config["ref"]["fasta"],
        ref_dict=config["ref"]["dict"],
        bam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_sample(w, file="bed_padded")
    output:
        "%s/calling/germline_haplotype_caller/{nsample}.vcf.gz" % R_FOLDER
    benchmark:
        "%s/calling/germline_haplotype_caller/{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/germline_haplotype_caller/{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        stand_min_conf=config["params"]["gatk"]["HaplotypeCaller"]["stand_min_conf"]
    threads: 4
    resources:
        queue="mediumq",
        mem_mb=30000,
        time_min=1200
    shell:
        """gatk --java-options {params.java} HaplotypeCaller \
            --reference {input.ref} \
            --native-pair-hmm-threads {threads} \
            --standard-min-confidence-threshold-for-calling {params.stand_min_conf} \
            --input {input.bam} \
            --output {output} 2> {log}"""


# Use hard-threshold heuristics to filter-out putative false-positives. Filters and thresholds differ slightly between
# snps and indels. Create a filtered VCF in which passing variants are annotated as PASS and failing variants are
# annotated with the name(s) of the filter(s) they failed.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variantsh
rule germline_filter_false_positives:
    input:
        "%s/calling/germline_haplotype_caller/{nsample}.vcf.gz" % R_FOLDER
    output:
        snp=temp("%s/calling/germline_filter_false_positives/{nsample}_snp.vcf.gz" % R_FOLDER),
        snp_tbi=temp("%s/calling/germline_filter_false_positives/{nsample}_snp.vcf.gz.tbi" % R_FOLDER),
        indel=temp("%s/calling/germline_filter_false_positives/{nsample}_indel.vcf.gz" % R_FOLDER),
        indel_tbi=temp("%s/calling/germline_filter_false_positives/{nsample}_indel.vcf.gz.tbi" % R_FOLDER),
        snp_filt=temp("%s/calling/germline_filter_false_positives/{nsample}_snp_filt.vcf.gz" % R_FOLDER),
        snp_filt_tbi=temp("%s/calling/germline_filter_false_positives/{nsample}_snp_filt.vcf.gz.tbi" % R_FOLDER),
        indel_filt=temp("%s/calling/germline_filter_false_positives/{nsample}_indel_filt.vcf.gz" % R_FOLDER),
        indel_filt_tbi=temp("%s/calling/germline_filter_false_positives/{nsample}_indel_filt.vcf.gz.tbi" % R_FOLDER),
        vcf_filt="%s/calling/germline_filter_false_positives/{nsample}.vcf.gz" % R_FOLDER,
        vcf_filt_tbi="%s/calling/germline_filter_false_positives/{nsample}.vcf.gz.tbi" % R_FOLDER
    benchmark:
        "%s/calling/germline_filter_false_positives/{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/germline_filter_false_positives/{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx10g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        snp_filters='-filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" ' +
            '-filter \"SOR > 3.0\" --filter-name \"SOR3\" -filter \"FS > 60.0\" --filter-name \"FS60\" ' +
            '-filter \"MQ < 40.0\" --filter-name \"MQ40\" -filter \"MQRankSum < -12.5\" --filter-name \"MQRankSum-12.5\" '
            '-filter \"ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSum-8\"',
        indel_filters=' -filter \"QD < 2.0\" --filter-name \"QD2\" -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" '+
            '-filter \"FS > 200.0\" --filter-name \"FS200\" -filter \"ReadPosRankSum < -20.0\" ' +
            '--filter-name \"ReadPosRankSum-20\"'
    threads: 1
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=60
    shell:
        "gatk --java-options {params.java} SelectVariants -V {input} -select-type SNP -O {output.snp} 2>> {log} &&"
        "gatk --java-options {params.java} SelectVariants -V {input} -select-type MIXED -select-type INDEL -O {output.indel} 2>> {log} &&"
        "gatk --java-options {params.java} VariantFiltration -V {output.snp} {params.snp_filters} -O {output.snp_filt} 2>> {log} &&"
        "gatk --java-options {params.java} VariantFiltration -V {output.indel} {params.indel_filters} -O {output.indel_filt} 2>> {log} &&"
        "gatk --java-options {params.java} MergeVcfs -I {output.indel_filt} -I {output.snp_filt} -O {output.vcf_filt} 2>> {log}"


####
#### Somatic mutations ####
####


# Apply gatk mutect2 to call somatic variants. The following set of rules implement the mutect2 Read Orientation
# Artifacts Workflow and Panel of Normal workflow
# Tumor-Normal mode
rule somatic_mutect2_tumor_normal:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples]),
        nsample = "|".join([re.escape(x) for x in nsamples])
    input:
        ref=config["ref"]["fasta"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        nbam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        nbai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_tumor_normal(w, file="bed_padded"),
        pon="%s/calling/pon_merge/pon.vcf.gz" % R_FOLDER,
        germline_resource=config["params"]["gatk"]["Mutect2"]["germline_resource"]
    output:
        vcf="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        idx="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER,
        sta="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.vcf.gz.stats" % R_FOLDER,
        f1r2="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.f1r2.tar.gz" % R_FOLDER
    benchmark:
        "%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 4
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=240
    shell:
        """gatk --java-options {params.java} Mutect2 \
            --reference {input.ref} \
            --intervals {input.bed} \
            --native-pair-hmm-threads {threads} \
            --germline-resource {input.germline_resource} \
            --input {input.tbam} \
            --input {input.nbam} \
            --tumor-sample {wildcards.tsample} \
            --normal-sample {wildcards.nsample} \
            --panel-of-normals {input.pon} \
            --dont-use-soft-clipped-bases true \
            --f1r2-tar-gz {output.f1r2} \
            --output {output.vcf} 2> {log}"""


# Apply gatk mutect2 to call somatic variants. The following set of rules implement the mutect2 Read Orientation
# Artifacts Workflow and Panel of Normal workflow
# Tumor-only mode
rule somatic_mutect2_tumor_only:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples])
    input:
        ref=config["ref"]["fasta"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_tumor_normal(w, file="bed_padded"),
        pon="%s/calling/pon_merge/pon.vcf.gz" % R_FOLDER,
        germline_resource=config["params"]["gatk"]["Mutect2"]["germline_resource"]
    output:
        vcf="%s/calling/somatic_mutect2/{tsample}_vs_NA.vcf.gz" % R_FOLDER,
        idx="%s/calling/somatic_mutect2/{tsample}_vs_NA.vcf.gz.tbi" % R_FOLDER,
        sta="%s/calling/somatic_mutect2/{tsample}_vs_NA.vcf.gz.stats" % R_FOLDER,
        f1r2="%s/calling/somatic_mutect2/{tsample}_vs_NA.f1r2.tar.gz" % R_FOLDER
    benchmark:
        "%s/calling/somatic_mutect2/{tsample}_vs_NA.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_mutect2/{tsample}_vs_NA.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 4
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=240
    shell:
        """gatk --java-options {params.java} Mutect2 \
            --reference {input.ref} \
            --intervals {input.bed} \
            --native-pair-hmm-threads {threads} \
            --germline-resource {input.germline_resource} \
            --input {input.tbam} \
            --tumor-sample {wildcards.tsample} \
            --panel-of-normals {input.pon} \
            --dont-use-soft-clipped-bases true \
            --f1r2-tar-gz {output.f1r2} \
            --output {output.vcf} 2> {log}"""


# Part of mutect2 Read Orientation Artifacts Workflow and Panel of Normal workflow
# Learn read orientation model for filtering variants produced by Mutect2 in a later rule.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
rule somatic_learn_read_orientation_model:
    input:
        "%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.f1r2.tar.gz" % R_FOLDER
    output:
        temp("%s/calling/somatic_learn_read_orientation_model/{tsample}_vs_{nsample}_read-orientation_model.tar.gz" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_learn_read_orientation_model/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_learn_read_orientation_model/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx15g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=20000,
        time_min=20
    shell:
        """gatk --java-options {params.java} LearnReadOrientationModel \
            --input {input} \
            --output {output} 2> {log}"""


# Part of mutect2 Read Orientation Artifacts Workflow and Panel of Normal workflow
# Get pileup summaries table that will serve as input to compute cross-sample contamination
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
rule somatic_get_pileup_summaries:
    input:
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        variant=config["params"]["gatk"]["GetPileupSummaries"]["variant"]
    output:
        temp("%s/calling/somatic_get_pileup_summaries/{tsample}_get_pileup_summaries.table" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_get_pileup_summaries/{tsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_get_pileup_summaries/{tsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx20g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=20000,
        time_min=20
    shell:
        """gatk --java-options {params.java} GetPileupSummaries \
            --input {input.tbam} \
            --variant {input.variant} \
            --intervals {input.variant} \
            --output {output} 2> {log}"""


# Part of mutect2 Read Orientation Artifacts Workflow and Panel of Normal workflow
# Estimate cross-sample contamination
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
rule somatic_calculate_contamination:
    input:
        "%s/calling/somatic_get_pileup_summaries/{tsample}_get_pileup_summaries.table" % R_FOLDER
    output:
        temp("%s/calling/somatic_calculate_contamination/{tsample}_calculate_contamination.table" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_calculate_contamination/{tsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_calculate_contamination/{tsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """gatk --java-options {params.java} CalculateContamination \
            --input {input} \
            --output {output} 2> {log}"""


# Part of mutect2 Read Orientation Artifacts Workflow and Panel of Normal workflow
# Estimate cross-sample contamination
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2
rule somatic_maf_filter_mutect_calls:
    input:
        ref=config["ref"]["fasta"],
        vcf="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        sta="%s/calling/somatic_mutect2/{tsample}_vs_{nsample}.vcf.gz.stats" % R_FOLDER,
        rom="%s/calling/somatic_learn_read_orientation_model/{tsample}_vs_{nsample}_read-orientation_model.tar.gz" % R_FOLDER,
        con="%s/calling/somatic_calculate_contamination/{tsample}_calculate_contamination.table" % R_FOLDER
    output:
        vcf="%s/calling/somatic_maf_filter_mutect_calls/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_maf_filter_mutect_calls/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER,
        sta="%s/calling/somatic_maf_filter_mutect_calls/{tsample}_vs_{nsample}.vcf.gz.filteringStats.tsv" % R_FOLDER
    benchmark:
        "%s/calling/somatic_maf_filter_mutect_calls/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_filter_mutect_calls/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=60
    shell:
        """gatk --java-options {params.java} FilterMutectCalls \
            --reference {input.ref} \
            --variant {input.vcf} \
            --orientation-bias-artifact-priors {input.rom} \
            --contamination-table {input.con} \
            --output {output.vcf} 2> {log}"""


# For TCGA, build vcf files with MC3 filters using the maf2vcf tool.
rule somatic_maf_convert_maf2vcf:
    input:
        ref=config["ref"]["fasta"],
        maf="%s/data/tcga_pancanatlas/mc3.v0.2.8.CONTROLLED.maf.gz" % R_FOLDER
    output:
        maf=temp("%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.maf" % R_FOLDER),
        tsv=temp("%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.pairs.tsv" % R_FOLDER),
        vcf=temp("%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.vcf" % R_FOLDER),
        gzvcf="%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        gztbi="%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER
    benchmark:
        "%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_convert_maf2vcf/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    params:
        awk=lambda w: "'NR==1; NR>1{if ($16==\"%s\" && $17==\"%s\") {print}}'" % (w.tsample, w.nsample),
        vcf2maf=config["params"]["vcf2maf"]["path"],
        dir="%s/calling/somatic_maf_convert_maf2vcf" % R_FOLDER
    threads: 4
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=60
    shell:
        "zcat {input.maf} | awk -F '\\t' {params.awk} > {output.maf} && "
        "python -u workflow/scripts/utils_prepare_convert_maf2vcf.py --input {output.maf} --output {output.maf} && "
        "perl {params.vcf2maf}/maf2vcf.pl --input-maf {output.maf} --output-dir {params.dir} --ref-fasta {input.ref} && "
        "bgzip -c {output.vcf} > {output.gzvcf} && "
        "tabix -p vcf {output.gzvcf}"


# Use hard-threshold heuristics to filter-out putative false-positives. Filters and thresholds differ slightly between
# snps and indels. Create a filtered VCF in which passing variants are annotated as PASS and failing variants are
# annotated with the name(s) of the filter(s) they failed.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variantsh
rule somatic_maf_filter_false_positives:
    input:
        get_input_somatic_maf_filter_false_positives
    output:
        vcf=temp("%s/calling/somatic_maf_filter_false_positives/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER),
        tbi=temp("%s/calling/somatic_maf_filter_false_positives/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_maf_filter_false_positives/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_filter_false_positives/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx10g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        filters='-filter \"vc.getGenotype(\'{nsample}\').getDP() < 10\" --filter-name \"LOW_COVERAGE_NORMAL\" ' +
            '-filter \"vc.getGenotype(\'{tsample}\').getDP() < 20\" --filter-name \"LOW_COVERAGE_TUMOR\" ' +
            '-filter \"(vc.getGenotype(\'{tsample}\').getAD().1)/(1.0*vc.getGenotype(\'{tsample}\').getDP()) < 0.05\" --filter-name \"LOW_VAF\"'
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=20
    shell:
        """gatk --java-options {params.java} VariantFiltration \
            --variant {input} \
            {params.filters} \
            --output {output.vcf} 2> {log}"""


# Add a filter for variants not in the intersection of all target files of the project.
rule somatic_maf_filter_outside_intersection:
    input:
        vcf="%s/calling/somatic_maf_filter_false_positives/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        bed=config["target_files"]["bed_padded"]["all_targets_intersect"]
    output:
        vcf=temp("%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER),
        tbi=temp("%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        dir="%s/calling/somatic_maf_filter_outside_intersection" % R_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=20
    shell:
        """
        bash workflow/scripts/03.1_maf_filter_outside_intersection.sh \
            -d {params.dir} \
            -v {input.vcf} \
            -t {wildcards.tsample} \
            -n {wildcards.nsample} \
            -b {input.bed} \
            -o {output.vcf} 2> {log}"""


# Save VCF after all filtering.
# Extract a tab-delimited format file from the VCF with minimal information about the filters applied on variants.
rule somatic_maf_filters_vcf_to_table:
    input:
        vcf="%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER
    output:
        tsv_tmp_1=temp("%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}_tmp_1.tsv" % R_FOLDER),
    benchmark:
        "%s/calling/somatic_maf_filters_vcf_to_table/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_filters_vcf_to_table/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    params:
        java="'-Xmx4g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        gatk --java-options {params.java} VariantsToTable \
            -V {input.vcf} \
            -F CHROM -F POS -F REF -F ALT -F TYPE -F FILTER \
            --show-filtered \
            -O {output.tsv_tmp_1} 2> {log}
        """

rule somatic_maf_filters_process_fields:
    input:
        vcf="%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_maf_filter_outside_intersection/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER,
        tsv_tmp_1="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}_tmp_1.tsv" % R_FOLDER
    output:
        vcf="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER,
        tsv_tmp_2=temp("%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}_tmp_2.tsv" % R_FOLDER),
        tsv="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER
    benchmark:
        "%s/calling/somatic_maf_filters_process_fields/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_filters_process_fields/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        cp {input.vcf} {output.vcf} && \
        cp {input.tbi} {output.tbi} && \
        python -u workflow/scripts/utils_add_sample_ids.py \
            --input {input.tsv_tmp_1} \
            --tsample {wildcards.tsample} \
            --tlabel Tumor_Sample \
            --nsample {wildcards.nsample} \
            --nlabel Normal_Sample \
            --output {output.tsv_tmp_2} 2> {log} && \
        python -u workflow/scripts/utils_remove_empty_fields.py  \
            --input {output.tsv_tmp_2} \
            --level 1 \
            --output {output.tsv} 2>> {log}
        """


# Extract only PASS variants for annotation afterwards.
rule somatic_maf_select_pass:
    input:
        vcf="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_maf_filters/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER
    output:
        "%s/calling/somatic_maf_select_pass/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER
    benchmark:
        "%s/calling/somatic_maf_select_pass/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_maf_select_pass/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    params:
        java="'-Xmx4g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=10
    shell:
        """
        gatk --java-options {params.java} SelectVariants \
            --variant {input.vcf} \
            --exclude-filtered \
            --output {output} 2> {log}
        """

####
#### Copy number variants ####
####

# SNP pileup table for CNV calling
rule somatic_cnv_get_snp_pileup_tumor_normal:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples]),
        nsample = "|".join([re.escape(x) for x in nsamples])
    log:
        "%s/calling/somatic_cnv_get_snp_pileup/{tsample}_vs_{nsample}.log" % L_FOLDER
    benchmark:
        "%s/calling/somatic_cnv_get_snp_pileup/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    conda:
        "../envs/main.yaml"
    input:
        vcf=config["params"]["gatk"]["BaseRecalibrator"]["known_sites"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        nbam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        nbai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER
    output:
        snp_pileup="%s/calling/somatic_snp_pileup/{tsample}_vs_{nsample}.csv.gz" % R_FOLDER,
        nbhd_snp="%s/calling/somatic_nbhd_snp/{tsample}_vs_{nsample}.tsv" % R_FOLDER
    threads: 10
    resources:
        queue="shortq",
        mem_mb=64000,
        time_min=60
    shell:
        """
        Rscript workflow/scripts/05.1_get_snp_pileup.R \
          -t {input.tbam} \
          -n {input.nbam} \
          -op {output.snp_pileup} \
          -on {output.nbhd_snp} \
          -vcf {input.vcf} \
          -N {threads} &> {log}
        """


# Call CNV using cnv-facets which is a wrapper around facets
# See https://github.com/dariober/cnv_facets
# Tumor-normal mode
rule somatic_cnv_facets_tumor_normal:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples]),
        nsample = "|".join([re.escape(x) for x in nsamples])
    input:
        snp_pileup="%s/calling/somatic_snp_pileup/{tsample}_vs_{nsample}.csv.gz" % R_FOLDER,
        nbhd_snp="%s/calling/somatic_nbhd_snp/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        bed=lambda w: get_target_file_tumor_normal(w, file="bed_padded")
    output:
        vcf="%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.vcf.gz.tbi" % R_FOLDER,
        png="%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.cnv.png" % R_FOLDER,
        spider="%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.spider.pdf" % R_FOLDER
    benchmark:
        "%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        code_dir=config["params"]["cnv"]["cnv_facets"]["code_dir"],
        diplogr=get_facets_diplogr,
        prefix="{tsample}_vs_{nsample}",
        nbhd_snp=lambda wildcards, input: pd.read_table(input.nbhd_snp)["nbhd_snp"].get(0),
        cval_pre=config["params"]["cnv"]["facets"]["cvals"]["pre"],
        cval_pro=config["params"]["cnv"]["facets"]["cvals"]["pro"],
        gbuild=config["params"]["cnv"]["facets"]["gbuild"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=64000,
        time_min=60
    shell:
        """
        Rscript {params.code_dir}/bin/cnv_facets.R \
            -p {input.snp_pileup} \
            -snp {params.nbhd_snp} \
            -o {params.prefix} \
            -N {threads} \
            -T {input.bed} \
            -r {params.diplogr} \
            -np \
            --gbuild {params.gbuild} \
            --cval {params.cval_pre} {params.cval_pro} &> {log} && \
        mv {wildcards.tsample}_vs_{wildcards.nsample}.vcf.gz {output.vcf} && \
        mv {wildcards.tsample}_vs_{wildcards.nsample}.vcf.gz.tbi {output.tbi} && \
        mv {wildcards.tsample}_vs_{wildcards.nsample}.cnv.png {output.png} && \
        mv {wildcards.tsample}_vs_{wildcards.nsample}.spider.pdf {output.spider}
	"""


# SNP pileup table for CNV calling
rule somatic_cnv_get_snp_pileup_tumor_unmatched:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples])
    log:
        "%s/somatic_cnv_get_snp_pileup/{tsample}_vs_%s.log" % (L_FOLDER, nsample_normal_pool)
    benchmark:
        "%s/somatic_cnv_get_snp_pileup/{tsample}_vs_%s.tsv" % (B_FOLDER, nsample_normal_pool)
    conda:
        "../envs/main.yaml"
    input:
        vcf=config["params"]["gatk"]["BaseRecalibrator"]["known_sites"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        nbam="%s/mapping/%s.nodup.recal.bam" % (R_FOLDER, nsample_normal_pool),
        nbai="%s/mapping/%s.nodup.recal.bam.bai" % (R_FOLDER, nsample_normal_pool)
    output:
        snp_pileup="%s/calling/somatic_snp_pileup/{tsample}_vs_%s.csv.gz" % (R_FOLDER, nsample_normal_pool),
        nbhd_snp="%s/calling/somatic_nbhd_snp/{tsample}_vs_%s.tsv" % (R_FOLDER, nsample_normal_pool)
    threads: 10
    resources:
        queue="shortq",
        mem_mb=64000,
        time_min=60
    shell:
        """
        Rscript workflow/scripts/05.1_get_snp_pileup.R \
          -t {input.tbam} \
          -n {input.nbam} \
          -op {output.snp_pileup} \
          -on {output.nbhd_snp} \
          -vcf {input.vcf} \
          -N {threads} &> {log}
        """


# Call CNV using cnv-facets which is a wrapper around facets
# See https://github.com/dariober/cnv_facets
# Tumor-unmatched mode
rule somatic_cnv_facets_tumor_unmatched:
    wildcard_constraints:
        tsample = "|".join([re.escape(x) for x in tsamples])
    input:
        vcf=config["params"]["gatk"]["BaseRecalibrator"]["known_sites"],
        snp_pileup="%s/calling/somatic_snp_pileup/{tsample}_vs_%s.csv.gz" % (R_FOLDER, nsample_normal_pool),
        nbhd_snp="%s/calling/somatic_nbhd_snp/{tsample}_vs_%s.tsv" % (R_FOLDER, nsample_normal_pool)
    output:
        vcf="%s/calling/somatic_cnv_facets/{tsample}_vs_NA.vcf.gz" % R_FOLDER,
        tbi="%s/calling/somatic_cnv_facets/{tsample}_vs_NA.vcf.gz.tbi" % R_FOLDER,
        png="%s/calling/somatic_cnv_facets/{tsample}_vs_NA.cnv.png" % R_FOLDER,
        cov="%s/calling/somatic_cnv_facets/{tsample}_vs_NA.cov.pdf" % R_FOLDER,
        spider="%s/calling/somatic_cnv_facets/{tsample}_vs_NA.spider.pdf" % R_FOLDER
    benchmark:
        "%s/calling/somatic_cnv_facets/{tsample}_vs_NA.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_cnv_facets/{tsample}_vs_NA.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        code_dir=config["params"]["cnv"]["cnv_facets"]["code_dir"],
        prefix="{tsample}_vs_NA",
        nbhd_snp=lambda wildcards, input: pd.read_table(input.nbhd_snp)["nbhd_snp"].get(0),
        cval_pre=config["params"]["cnv"]["facets"]["cvals"]["pre"],
        cval_pro=config["params"]["cnv"]["facets"]["cvals"]["pro"],
        gbuild=config["params"]["cnv"]["facets"]["gbuild"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=60000,
        time_min=60
    shell:
        """
        Rscript {params.code_dir}/bin/cnv_facets.R \
            -vcf {input.vcf} \
            -p {input.snp_pileup} \
            -snp {params.nbhd_snp} \
            -u \
            -o {params.prefix} \
            -N {threads} \
            --gbuild {params.gbuild} \
            --cval {params.cval_pre} {params.cval_pro} &> {log} && \
        mv {wildcards.tsample}_vs_NA.vcf.gz {output.vcf} && \
        mv {wildcards.tsample}_vs_NA.vcf.gz.tbi {output.tbi} && \
        mv {wildcards.tsample}_vs_NA.cnv.png {output.png} && \
        mv {wildcards.tsample}_vs_NA.cov.pdf {output.cov} && \
        mv {wildcards.tsample}_vs_NA.spider.pdf {output.spider}
        """


# Process VCF file
rule somatic_cnv_process_vcf:
    input:
        vcf="%s/calling/somatic_cnv_facets/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        rules_arm=config["params"]["cnv"]["chr_arm_rules"],
        rules_cat=config["params"]["cnv"]["cna_cat_rules"],
        env="%s/setup_r.done" % L_FOLDER
    output:
        arm="%s/calling/somatic_cnv_chr_arm/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        sum="%s/calling/somatic_cnv_sum/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        tab="%s/calling/somatic_cnv_table/{tsample}_vs_{nsample}.tsv" % R_FOLDER
    benchmark:
        "%s/calling/somatic_cnv_process_vcf/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_cnv_process_vcf/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/r.yaml"
    threads: 1
    params:
        gender = lambda w: get_column_table_sample(w, "Gender"),
        threshold=config["params"]["cnv"]["calls_threshold"]
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=20
    shell:
        """
        Rscript workflow/scripts/05.2_cnv_process_vcf.R \
            --input_vcf {input.vcf} \
            --gender {params.gender} \
            --rules_arm {input.rules_arm} \
            --rules_cat {input.rules_cat} \
            --threshold {params.threshold} \
            --output_arm {output.arm} \
            --output_sum {output.sum} \
            --output_tab {output.tab} \
            --log {log}
        """


# Convert table with cnv at segments to cnv at genes using bedtools
rule somatic_cnv_gene_calls:
    input:
        tab="%s/calling/somatic_cnv_table/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        bed=config["params"]["cnv"]["bed"]
    output:
        "%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER
    benchmark:
        "%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    params:
        threshold=config["params"]["cnv"]["calls_threshold"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=20
    shell:
        """
        python -u workflow/scripts/05.3_cnv_gene_calls.py \
            --input_tab {input.tab} \
            --input_bed {input.bed} \
            --threshold {params.threshold} \
            --output {output} &> {log}
        """


# Make a table of filter cnv calls per gene
rule somatic_cnv_gene_calls_filtered:
    input:
        "%s/calling/somatic_cnv_gene_calls/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER
    output:
        temp("%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER)
    benchmark:
        "%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.log" % L_FOLDER
    threads: 1
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=20
    shell:
        """
        zcat {input} | grep "PASS\|Tumor_Sample_Barcode" | gzip > {output} 2> {log}
        """


####
#### Microsatellite instability ####
####


# Call MSI using MANTIS
# See https://github.com/OSU-SRLab/MANTIS
rule somatic_msi_mantis:
    input:
        ref=config["ref"]["fasta"],
        code="%s/mantis.py" % config["params"]["mantis"]["code_dir"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        nbam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        nbai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER,
        ms="%s/%s_microsatellites_%s" %
                (config["params"]["mantis"]["ms_dir"],
                 config["params"]["mantis"]["buildver"],
                 os.path.basename(config["params"]["mantis"]["ms_bed"]))
    output:
        msi="%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        msi_status="%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.tsv.status" % R_FOLDER,
        kmer_counts="%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.kmer_counts.tsv" % R_FOLDER,
        kmer_counts_filtered="%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.kmer_counts_filtered.tsv" % R_FOLDER
    benchmark:
        "%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/somatic_msi_mantis/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    params:
        mrq=20,
        mlq=25,
        mlc=20,
        mrr=1
    threads: 4
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=30
    shell:
        """
        {input.code} \
            --genome {input.ref} \
            --bedfile {input.ms} \
            --tumor {input.tbam} \
            --normal {input.nbam} \
            -mrq {params.mrq} \
            -mlq {params.mlq} \
            -mlc {params.mlc} \
            -mrr {params.mrr} \
            --threads {threads} \
            --output {output.msi} &> {log}
        """
