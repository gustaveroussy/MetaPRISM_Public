####
#### TCGA specific ####
####

# Download TCGA data files from manifests using gdc-client.
rule download_gdc_tcga_pancanatlas:
    input:
        client=config["download_gdc"]["client"],
        manifest=config["download_gdc"]["manifests"]["tcga_pancanatlas"],
        token=config["download_gdc"]["token"]
    output:
        expand("%s/data/tcga_pancanatlas/{id}" % R_FOLDER, id=ids_gdc_tcga_pancanatlas)
    benchmark:
        "%s/download/download_gdc_tcga_pancanatlas.tsv" % B_FOLDER
    log:
        "%s/download/download_gdc_tcga_pancanatlas.log" % L_FOLDER
    params:
        dir="%s/data/tcga_pancanatlas" % R_FOLDER
    threads: 16
    resources:
        queue="shortq",
        mem_mb=24000,
        time_min=240
    shell:
        """
        ./{input.client} download \
            --manifest {input.manifest} \
            --token-file {input.token} \
            --dir {params.dir} \
            --n-processes {threads} \
            --log-file {log}
        """


# Download BAM files from manifests using gdc-client.
rule download_gdc_bam:
    log:
        "%s/fastq/gdc_get_bam_{bam_name}.log" % L_FOLDER
    benchmark:
        "%s/fastq/gdc_get_bam_{bam_name}.tsv" % B_FOLDER
    input:
        table = config["samples"],
        gdc_client = config["params"]["gdc"]["client"],
        gdc_token = config["params"]["gdc"]["token"]
    params:
        dir = "%s/data/bam" % R_FOLDER,
        bam_id = get_bam_id
    output:
        "%s/data/bam/{bam_name}" % R_FOLDER
    resources:
        load=40,
        queue="mediumq",
        mem_mb=1000,
        time_min=960
    threads: 4
    shell:
        """
        bash workflow/scripts/01.3_gdc_get_bam.sh \
            -d {params.dir} \
            -i {params.bam_id} \
            -n {wildcards.bam_name} \
            -c {input.gdc_client} \
            -t {input.gdc_token} \
            -a {threads} &> {log}
        """


# Convert BAM to FASTQs.
rule bam_to_fastq:
    log:
        "%s/fastq/bam_to_fastq_{sample}.log" % L_FOLDER
    benchmark:
        "%s/fastq/bam_to_fastq_{sample}.tsv" % B_FOLDER
    input:
        lambda w: "%s/data/bam/%s" % (R_FOLDER, get_column_table_sample(w, "BAM_Name"))
    output:
        r1=temp("%s/data/fastq/{sample}_R1.fastq") % R_FOLDER,
        r2=temp("%s/data/fastq/{sample}_R2.fastq") % R_FOLDER
    conda:
        "../envs/main.yaml"
    threads: 1
    resources:
        queue="mediumq",
        mem_mb=1000,
        time_min=360
    shell:
        """
        bash workflow/scripts/01.4_bam_to_fastq.sh \
            -b {input} \
            -x {output.r1} \
            -y {output.r2} &> {log}
        """


# Compress FASTQs.
rule fastq_compress:
    log:
        "%s/fastq/fastq_compress_{sample}_{i}.log" % L_FOLDER
    benchmark:
        "%s/fastq/fastq_compress_{sample}_{i}.tsv" % B_FOLDER
    input:
        "%s/data/fastq/{sample}_R{i}.fastq" % R_FOLDER,
    output:
        "%s/data/fastq/{sample}_R{i}.fastq.gz" % R_FOLDER,
    threads: 1
    resources:
        queue="shortq",
        mem_mb=1000,
        time_min=120
    shell:
        """
        gzip {input}
        """
