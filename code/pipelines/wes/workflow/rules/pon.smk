# Apply gatk mutect2 to call variants from the normal samples in tumor-mode. The variants called 
# will serve as a basis for creating the panel-of-normals
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA
rule pon_mutect2:
    input:
        ref=config["ref"]["fasta"],
        nbam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        nbai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_sample(w, file="bed_padded"),
        germline_resource=config["params"]["gatk"]["Mutect2"]["germline_resource"]
    output:
        vcf="%s/calling/pon_mutect2/{nsample}.vcf.gz" % R_FOLDER,
        idx="%s/calling/pon_mutect2/{nsample}.vcf.gz.tbi" % R_FOLDER,
        sta="%s/calling/pon_mutect2/{nsample}.vcf.gz.stats" % R_FOLDER
    benchmark:
        "%s/calling/pon_mutect2/{nsample}.tsv" % B_FOLDER
    log:
        "%s/calling/pon_mutect2/{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx30g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 4
    resources:
        queue="shortq",
        mem_mb=20000,
        time_min=240
    shell:
        """gatk --java-options {params.java} Mutect2 \
            --reference {input.ref} \
            --intervals {input.bed} \
            --native-pair-hmm-threads {threads} \
            --germline-resource {input.germline_resource} \
            --max-mnp-distance 0 \
            --input {input.nbam} \
            --tumor-sample {wildcards.nsample} \
            --output {output.vcf} 2> {log}"""


# Apply gatk GenomicsDBImport before creating the panel-of-normals.
# See https://javadoc.io/static/org.broadinstitute/gatk/4.1.4.0/org/broadinstitute/hellbender/tools/walkers/mutect/CreateSomaticPanelOfNormals.html
rule pon_genomicsdb:
    input:
        ref=config["ref"]["fasta"],
        vcfs=get_input_vcfs_pon_genomicsdb,
        bed=get_input_bed_pon_genomicsdb
    output:
        temp(directory("%s/calling/pon_genomicsdb/{target}_{interval}" % R_FOLDER))
    benchmark:
        "%s/calling/pon_genomicsdb/{target}_{interval}.tsv" % B_FOLDER
    log:
        "%s/calling/pon_genomicsdb/{target}_{interval}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx25g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp -DGATK_STACKTRACE_ON_USER_EXCEPTION=true'",
        vcfs=lambda w: "-V " + " -V ".join(get_input_vcfs_pon_genomicsdb(w))
    threads: 1
    resources:
        queue="shortq",
        mem_mb=25000,
        time_min=120
    shell:
        """gatk --java-options {params.java} GenomicsDBImport \
            --reference {input.ref} \
            --intervals {input.bed} \
            {params.vcfs} \
            --genomicsdb-workspace-path {output} \
            --overwrite-existing-genomicsdb-workspace 2> {log}"""


# Apply gatk CreateSomaticPanelOfNormals to create the panel-of-normals vcf
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA
rule pon_create:
    input:
        ref=config["ref"]["fasta"],
        pon_db="%s/calling/pon_genomicsdb/{target}_{interval}" % R_FOLDER
    output:
        vcf=temp("%s/calling/pon_create/pon_{target}_{interval}.vcf.gz" % R_FOLDER),
        tbi=temp("%s/calling/pon_create/pon_{target}_{interval}.vcf.gz.tbi" % R_FOLDER)
    benchmark:
        "%s/calling/pon_create/{target}_{interval}.tsv" % B_FOLDER
    log:
        "%s/calling/pon_create/{target}_{interval}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx50g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        min_sample_count=config["params"]["gatk"]["CreateSomaticPanelOfNormals"]["min_sample_count"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=40000,
        time_min=45
    shell:
        """gatk --java-options {params.java} CreateSomaticPanelOfNormals \
            --reference {input.ref} \
            --variant gendb://{input.pon_db} \
            --min-sample-count {params.min_sample_count} \
            --output {output.vcf} 2> {log}
            rm -r {input.pon_db}
        """

# Sort VCFs created in CreateSomaticPanelOfNormals before merging them. This is required
# as MergeVcfs expects sorted vcfs in input.
rule pon_sort:
    input:
        vcf="%s/calling/pon_create/pon_{target}_{interval}.vcf.gz" % R_FOLDER,
        tbi="%s/calling/pon_create/pon_{target}_{interval}.vcf.gz.tbi" % R_FOLDER
    output:
        vcf=temp("%s/calling/pon_sort/pon_sorted_{target}_{interval}.vcf.gz" % R_FOLDER),
        tbi=temp("%s/calling/pon_sort/pon_sorted_{target}_{interval}.vcf.gz.tbi" % R_FOLDER)
    benchmark:
        "%s/calling/pon_sort/{target}_{interval}.tsv" % B_FOLDER
    log:
        "%s/calling/pon_sort/{target}_{interval}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx10g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=15
    shell:
        """gatk --java-options {params.java} SortVcf \
            --INPUT {input.vcf} \
            --OUTPUT {output.vcf} 2> {log}"""


# Apply gatk MergeVcfs on the vcfs created for the pon.
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037226612-MergeVcfs-Picard-
rule pon_merge:
    input:
        vcf=expand("%s/calling/pon_sort/pon_sorted_{target}_{interval}.vcf.gz" % R_FOLDER, get_allowed_pairs_target_interval(),
             target=targets, interval=intervals),
        tbi=expand("%s/calling/pon_sort/pon_sorted_{target}_{interval}.vcf.gz.tbi" % R_FOLDER, get_allowed_pairs_target_interval(),
             target=targets, interval=intervals)
    output:
        "%s/calling/pon_merge/pon.vcf.gz" % R_FOLDER
    benchmark:
        "%s/calling/pon_merge/pon.tsv" % B_FOLDER
    log:
        "%s/calling/pon_merge/pon.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx10g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        input_formatted=lambda wildcards, input: "-I " + " -I ".join(input.vcf)
    threads: 1
    resources:
        queue="shortq",
        mem_mb=10000,
        time_min=15
    shell:
        """gatk --java-options {params.java} MergeVcfs \
            {params.input_formatted} \
            --OUTPUT {output} 2> {log}"""
