# Align the fastq  file to the reference fastq file and sort with samtools to create the bam file
# See http://bio-bwa.sourceforge.net/bwa.shtml
rule bwa_mem:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    input:
        config["ref"]["fasta"],
        '%s/data/fastp/{sample}_R1.fastq.gz' % R_FOLDER,
        '%s/data/fastp/{sample}_R2.fastq.gz' % R_FOLDER
    output:
        temp("%s/mapping/{sample}.bam" % R_FOLDER)
    benchmark:
        "%s/mapping/bwa_mem/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bwa_mem/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        extra="-M -R '@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:truseq'"
    threads: 8
    resources:
        queue = "shortq",
        mem_mb = 50000,
        time_min = 360
    shell:
        """(bwa mem -t {threads} {params.extra} {input} \
            | samtools view -@ {threads} -b - \
            | samtools sort -@ {threads} - -o {output}) 2> {log}"""


# Remove duplicated reads from the bam file with MarkDuplicates (picard)
# See https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard- 
rule mark_duplicates:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    input:
        "%s/mapping/{sample}.bam" % R_FOLDER
    output:
        bam=temp("%s/mapping/{sample}.nodup.bam" % R_FOLDER),
        metrics="%s/qc/mark_duplicates/{sample}.metrics.txt" % R_FOLDER
    benchmark:
        "%s/mapping/mark_duplicates/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/mark_duplicates/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    params:
        java="'-Xmx120g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
        remove_duplicates=config["params"]["picard"]["MarkDuplicates"]["remove_duplicates"]
    resources:
        queue="shortq",
        mem_mb=120000,
        time_min=120
    shell:
        """gatk --java-options {params.java} MarkDuplicates \
            --INPUT {input} \
            --REMOVE_DUPLICATES {params.remove_duplicates} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics} 2> {log}"""


# Generate bam index with samtools before recalibration
# See http://www.htslib.org/doc/samtools-index.html
rule index_bam_before_recal:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    input:
        "%s/mapping/{sample}.nodup.bam" % R_FOLDER,
    output:
        temp("%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER)
    benchmark:
        "%s/mapping/index_bam_before_recal/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/index_bam_before_recal/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=300,
        time_min=20
    shell:
        """samtools index {input} 2> {log}"""


# Recalibrate base qualities using gatk BaseRecalibrator - first pass
rule bqsr_1:
    input:
        ref=config["ref"]["fasta"],
        ref_dict=config["ref"]["dict"],
        known_sites=config["params"]["gatk"]["BaseRecalibrator"]["known_sites"],
        bam="%s/mapping/{sample}.nodup.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_sample(w, file="bed_padded")
    output:
        temp("%s/mapping/{sample}_bqsr_1.table" % R_FOLDER)
    benchmark:
        "%s/mapping/bqsr_1/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/bqsr_1/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=50000,
        time_min=120
    shell:
        """gatk --java-options {params.java} BaseRecalibrator \
            --reference {input.ref} \
            --known-sites {input.known_sites} \
            --intervals {input.bed} \
            --input {input.bam} \
            --output {output} 2> {log}"""


# Recalibrate base qualities using gatk ApplyBQSR - first pass
rule apply_bqsr_1:
    input:
        bam="%s/mapping/{sample}.nodup.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.bam.bai" % R_FOLDER,
        ref=config["ref"]["fasta"],
        bqsr_recal_file="%s/mapping/{sample}_bqsr_1.table" % R_FOLDER
    output:
       temp("%s/mapping/{sample}.nodup.recal.beforeReformat.bam" % R_FOLDER)
    benchmark:
        "%s/mapping/apply_bqsr_1/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/apply_bqsr_1/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=15000,
        time_min=180
    shell:
        """gatk --java-options {params.java} ApplyBQSR \
            --create-output-bam-index false \
            --reference {input.ref} \
            --bqsr-recal-file {input.bqsr_recal_file} \
            --input {input.bam} \
            --output {output} 2> {log}"""


# Reformat bam files after base quality score recalibration with samtools in order to save space
rule reformat_bam_after_bqsr_1:
    input:
        "%s/mapping/{sample}.nodup.recal.beforeReformat.bam" % R_FOLDER
    output:
        "%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER
    benchmark:
        "%s/mapping/reformat_bam_after_bqsr_1/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/reformat_bam_after_bqsr_1/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 4
    resources:
        queue="shortq",
        mem_mb=10000,
        time_min=60
    shell:
        """samtools view \
            -@ {threads} \
            -b \
            -h \
            -o {output} {input} 2> {log}"""


# Generate bam index with samtools after recalibration
# See http://www.htslib.org/doc/samtools-index.html
rule index_bam_after_recal:
    wildcard_constraints:
        sample = "|".join([re.escape(x) for x in samples])
    input:
        "%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
    output:
        "%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
    benchmark:
        "%s/mapping/index_bam_after_recal/{sample}.tsv" % B_FOLDER
    log:
        "%s/mapping/index_bam_after_recal/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=300,
        time_min=20
    shell:
        """samtools index {input} 2> {log}"""


# Merge multiple normal bam files into one in order to use for unmatched normal
# CNV calling.
rule normal_pool_bam:
    input:
        expand("%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER, sample=nsamples_normal_pool)
    output:
        bam="%s/mapping/%s.nodup.recal.bam" % (R_FOLDER, nsample_normal_pool),
        idx="%s/mapping/%s.nodup.recal.bam.bai" % (R_FOLDER, nsample_normal_pool)
    benchmark:
        "%s/mapping/normal_pool_bam/%s.tsv" % (B_FOLDER, nsample_normal_pool)
    log:
        "%s/mapping/normal_pool_bam/%s.log" % (B_FOLDER, nsample_normal_pool)
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="mediumq",
        mem_mb=40000,
        time_min=240
    shell:
        """
        samtools merge {output.bam} {input} && \
        samtools index {output.bam} 2> {log}
        """
