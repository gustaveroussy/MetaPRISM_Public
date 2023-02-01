# A first run of fastqc on fastqs from IRODS
# See https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_1:
    input:
        '%s/data/fastq/{fastq}.fastq.gz' % R_FOLDER
    output:
        html='%s/qc/fastqc_1/{fastq}_fastqc.html' % R_FOLDER,
        zip='%s/qc/fastqc_1/{fastq}_fastqc.zip' % R_FOLDER
    benchmark:
        "%s/qc/fastqc_1/{fastq}.tsv" % B_FOLDER
    log:
        "%s/qc/fastqc_1/{fastq}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        adapters=config["params"]["fastqc"]["adapters"],
        out="%s/qc/fastqc_1" % R_FOLDER
    threads : 1
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=90
    shell:
        'fastqc -a {params.adapters} -o {params.out} {input} 2> {log}'


# Clean fastq files with fastp
# See https://github.com/OpenGene/fastp
rule fastp:
    input:
        unpack(get_fastqs_local)
    output:
        r1=temp('%s/data/fastp/{sample}_R1.fastq.gz' % R_FOLDER),
        r2=temp('%s/data/fastp/{sample}_R2.fastq.gz' % R_FOLDER),
        html='%s/qc/fastp/{sample}_fastp_report.html' % R_FOLDER,
        json='%s/qc/fastp/{sample}_fastp_report.json' % R_FOLDER
    benchmark:
        "%s/qc/fastp/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/fastp/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        adapters=config["params"]["fastp"]["adapters"],
        length_required=config["params"]["fastp"]["length_required"]
    threads : 12
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=240
    shell:
        """fastp --thread {threads} \
            --dont_overwrite \
            -i {input.r1} \
            -o {output.r1} \
            -I {input.r2} \
            -O {output.r2} \
            --compression 6 \
            --adapter_fasta {params.adapters} \
            --trim_poly_g \
            --trim_poly_x \
            --length_required {params.length_required} \
            --overrepresentation_analysis \
            --html {output.html} \
            --json {output.json} 2> {log}"""



# A second run of fastqc on fastqs processed by fastp
# See https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_2:
    input:
        '%s/data/fastp/{fastq}.fastq.gz' % R_FOLDER
    output:
        html='%s/qc/fastqc_2/{fastq}_fastqc.html' % R_FOLDER,
        zip='%s/qc/fastqc_2/{fastq}_fastqc.zip' % R_FOLDER
    benchmark:
        "%s/qc/fastqc_2/{fastq}.tsv" % B_FOLDER
    log:
        "%s/qc/fastqc_2/{fastq}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        adapters=config["params"]["fastqc"]["adapters"],
        out="%s/qc/fastqc_2" % R_FOLDER
    threads : 1
    resources:
        queue="shortq",
        mem_mb=5000,
        time_min=90
    shell:
        'fastqc -a {params.adapters} -o {params.out} {input} 2> {log}'


# Quality-control of mapping with gatk CollectHsMetrics (picard)
# See https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectHsMetrics-Picard-
rule collect_hs_metrics:
    input:
        bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER,
        intervals=lambda w: get_target_file_sample(w, file="intervals")
    output:
        "%s/qc/collect_hs_metrics/{sample}_hs_metrics.tsv" % R_FOLDER
    benchmark:
        "%s/qc/collect_hs_metrics/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/collect_hs_metrics/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        java="'-Xmx40g -Djava.io.tmpdir=/mnt/beegfs/scratch/tmp'",
    threads : 1
    resources:
        queue="shortq",
        mem_mb=30000,
        time_min=90
    shell:
        """gatk --java-options {params.java} CollectHsMetrics \
            --TI {input.intervals} \
            --BI {input.intervals} \
            --I {input.bam} \
            --O {output} 2> {log}"""


# Quality-control of mapping with samtools stats
# See http://www.htslib.org/doc/samtools-stats.html
rule samtools_stats:
    input:
        bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
    output:
        "%s/qc/samtools_stats/{sample}_stats.tsv" % R_FOLDER
    benchmark:
        "%s/qc/samtools_stats/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/samtools_stats/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads : 1
    resources:
        queue="shortq",
        mem_mb=3000,
        time_min=60
    shell:
        """samtools stats {input.bam} > {output} 2> {log}"""


# Quality-control of mapping with samtools flagstat
# See http://www.htslib.org/doc/samtools-flagstat.html
rule samtools_flagstat:
    input:
        bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER
    output:
        "%s/qc/samtools_flagstat/{sample}_flagstat.tsv" % R_FOLDER
    benchmark:
        "%s/qc/samtools_flagstat/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/samtools_flagstat/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    threads : 1
    resources:
        queue="shortq",
        mem_mb=3000,
        time_min=60
    shell:
        """samtools flagstat {input.bam} > {output} 2> {log}"""


# Quality-control of mapping with mosdepth
# See https://github.com/brentp/mosdepth
rule mosdepth:
    input:
        bam="%s/mapping/{sample}.nodup.recal.bam" % R_FOLDER,
        bai="%s/mapping/{sample}.nodup.recal.bam.bai" % R_FOLDER,
        bed=lambda w: get_target_file_sample(w, file="bed")
    output:
        "%s/qc/mosdepth/{sample}.mosdepth.global.dist.txt" % R_FOLDER
    benchmark:
        "%s/qc/mosdepth/{sample}.tsv" % B_FOLDER
    log:
        "%s/qc/mosdepth/{sample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        prefix="%s/qc/mosdepth/{sample}" % R_FOLDER
    threads : 4
    resources:
        queue="shortq",
        mem_mb=3000,
        time_min=60
    shell:
        """mosdepth --no-per-base \
            --threads {threads} \
            --fast-mode \
            --by {input.bed} \
            {params.prefix} {input.bam} 2> {log}"""


# Check that aligned sequencing files are derived from the same individual using NGSCheckMate.
rule ngs_checkmate_fastq:
    input:
        ref=config["ref"]["fasta"],
        snp_bed=config["params"]["ngs_checkmate"]["snp_bed"],
        tbam="%s/mapping/{tsample}.nodup.recal.bam" % R_FOLDER,
        tbai="%s/mapping/{tsample}.nodup.recal.bam.bai" % R_FOLDER,
        nbam="%s/mapping/{nsample}.nodup.recal.bam" % R_FOLDER,
        nbai="%s/mapping/{nsample}.nodup.recal.bam.bai" % R_FOLDER
    output:
        directory("%s/qc/ngs_checkmate/{tsample}_vs_{nsample}" % R_FOLDER)
    benchmark:
        "%s/qc/ngs_checkmate/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/qc/ngs_checkmate/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        code_dir=config["params"]["ngs_checkmate"]["code_dir"],
        prefix="{tsample}_vs_{nsample}"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=60
    shell:
        """
        export NCM_HOME={params.code_dir}
        echo {input.tbam} > {params.prefix}_bam_file_list.txt
        echo {input.nbam} >> {params.prefix}_bam_file_list.txt
        python -u {params.code_dir}/ncm.py \
            -B \
            -l {params.prefix}_bam_file_list.txt \
            -O {output} \
            -N {params.prefix} \
            -bed {input.snp_bed} &> {log}
        rm r_script.r.Rout
        rm {params.prefix}_bam_file_list.txt
        """
