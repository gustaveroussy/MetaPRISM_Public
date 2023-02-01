####
#### Germline mutations ####
####

# Run annovar on germline variants
# See https://annovar.openbioinformatics.org/en/latest/user-guide/startup/
rule germline_annovar:
    input:
        annovar=expand("%s/humandb/hg19_{resource}.txt" % config["params"]["annovar"]["path"],
            resource=config["params"]["annovar"]["protocol"]),
        vcf="%s/calling/germline_filter_false_positives/{nsample}.vcf.gz" % R_FOLDER
    output:
        "%s/annotation/germline_annovar/{nsample}.avinput" % R_FOLDER,
        "%s/annotation/germline_annovar/{nsample}.hg19_multianno.txt" % R_FOLDER,
        "%s/annotation/germline_annovar/{nsample}.hg19_multianno.vcf" % R_FOLDER
    benchmark:
        "%s/annotation/germline_annovar/{nsample}.tsv" % B_FOLDER
    log:
        "%s/annotation/germline_annovar/{nsample}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        annovar_dir=config["params"]["annovar"]["path"],
        output_dir="%s/annotation/annovar"  % R_FOLDER,
        buildver=config["params"]["annovar"]["buildver"],
        out="%s/annotation/germline_annovar/{nsample}" % R_FOLDER,
        protocol=",".join(config["params"]["annovar"]["protocol"]),
        operation=get_annovar_operation()
    threads: 4
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=240
    shell:
        """
        mkdir -p {params.output_dir}
        {params.annovar_dir}/table_annovar.pl {input.vcf} {params.annovar_dir}/humandb/ \
            -buildver {params.buildver} \
            --thread {threads} \
            --maxgenethread {threads} \
            -out {params.out} \
            -remove \
            -protocol {params.protocol} \
            -operation {params.operation} \
            -nastring . \
            -vcfinput 2> {log}"""


####
#### Somatic mutations ####
####

# Run vep on somatic variants
# See https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html#basic

# VEP with VCF output
rule somatic_vep_vcf:
    input:
        vcf="%s/calling/somatic_maf_select_pass/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
    output:
        vcf=temp("%s/annotation/somatic_vep/{tsample}_vs_{nsample}.vcf" % R_FOLDER),
    benchmark:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vep_vcf.tsv" % B_FOLDER
    log:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vep_vcf.log" % L_FOLDER
    conda:
        "../envs/perl.yaml"
    params:
        species=config["ref"]["species"],
        vep_dir=config["params"]["vep"]["path"],
        vep_data=config["params"]["vep"]["cache"],
    threads: 4
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=120
    shell:
        """
        {params.vep_dir}/vep  --input_file {input.vcf} \
            --dir {params.vep_data} \
            --cache \
            --canonical \
            --format vcf \
            --offline \
            --vcf \
            --no_progress \
            --no_stats \
            --af_gnomad \
            --appris \
            --sift b \
            --polyphen b \
            --biotype \
            --buffer_size 5000 \
            --hgvs \
            --species {params.species} \
            --symbol \
            --transcript_version \
            --output_file {output} \
            --warning_file {log}
        """


# VEP with Tab output
rule somatic_vep_tab:
    input:
        vcf="%s/calling/somatic_maf_select_pass/{tsample}_vs_{nsample}.vcf.gz" % R_FOLDER,
        cadd=config["params"]["vep"]["plugins_data"]["CADD"],
        dbnsfp=config["params"]["vep"]["plugins_data"]["dbNSFP"],
    output:
        vcf=temp("%s/annotation/somatic_vep/{tsample}_vs_{nsample}.tsv" % R_FOLDER),
    benchmark:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vep_tab.tsv" % B_FOLDER
    log:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vep_tab.log" % L_FOLDER
    conda:
        "../envs/perl.yaml"
    params:
        species=config["ref"]["species"],
        vep_dir=config["params"]["vep"]["path"],
        vep_data=config["params"]["vep"]["cache"],
        dbnsfp=",".join(config["params"]["vep"]["dbnsfp"])
    threads: 4
    resources:
        queue="shortq",
        mem_mb=10000,
        time_min=120
    shell:
        """
        {params.vep_dir}/vep  --input_file {input.vcf} \
            --dir {params.vep_data} \
            --cache \
            --canonical \
            --format vcf \
            --offline \
            --tab \
            --no_progress \
            --no_stats \
            --max_af \
            --af_gnomad \
            --sift b \
            --polyphen b \
            --appris \
            --biotype \
            --buffer_size 5000 \
            --hgvs \
            --mane \
            --species {params.species} \
            --symbol \
            --numbers \
            --transcript_version \
            --plugin CADD,{input.cadd[0]},{input.cadd[1]} \
            --plugin dbNSFP,{input.dbnsfp},{params.dbnsfp} \
            --output_file {output} \
            --warning_file {log}
        """


# Convert VCF annotated by VEP to MAF format
rule somatic_vep_vcf2maf:
    input:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}.vcf" % R_FOLDER
    output:
        temp("%s/annotation/somatic_vep/{tsample}_vs_{nsample}.maf" % R_FOLDER)
    benchmark:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vcf2maf.tsv" % B_FOLDER
    log:
        "%s/annotation/somatic_vep/{tsample}_vs_{nsample}_vcf2maf.log" % L_FOLDER
    threads: 1
    conda:
        "../envs/perl.yaml"
    params:
        path=config["params"]["vcf2maf"]["path"],
        fasta="%s/%s" % (config["params"]["vep"]["cache"],config["params"]["vep"]["fasta"]),
        dbnsfp=",".join(config["params"]["vep"]["dbnsfp"])
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        perl {params.path}/vcf2maf.pl \
            --inhibit-vep \
            --input-vcf {input} \
            --output {output} \
            --tumor-id {wildcards.tsample} \
            --normal-id {wildcards.nsample} \
            --ref-fasta {params.fasta} \
            --ncbi-build hg19 \
            --verbose 2> {log}
        """


# Merge the 2 tsv tables produced by VEP and by vcf2maf to produce final maf.
rule somatic_maf:
    input:
        vep="%s/annotation/somatic_vep/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        maf="%s/annotation/somatic_vep/{tsample}_vs_{nsample}.maf" % R_FOLDER
    output:
        "%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER
    benchmark:
        "%s/annotation/somatic_maf/{tsample}_vs_{nsample}.tsv" % B_FOLDER
    log:
        "%s/annotation/somatic_maf/{tsample}_vs_{nsample}.log" % L_FOLDER
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=20
    shell:
        """
        python -u workflow/scripts/03.2_maf_merge_vep_and_maf.py \
            --vep_table {input.vep} \
            --maf_table {input.maf} \
            --keep_vep_header \
            --output {output} &> {log}
        """


# Annnotate mutations using (in-house) civic annotator
if config["params"]["civic"]["run_per_sample"]["maf"]:
    # prepare a table for each pair tsample_vs_nsample
    rule somatic_maf_civic:
        input:
            table_alt="%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["civic"]["gene_list"],
            civic=config["params"]["civic"]["evidences"],
            rules=config["params"]["civic"]["rules_clean"]
        output:
            table_pre=temp("%s/annotation/somatic_maf_civic/{tsample}_vs_{nsample}_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_maf_civic/{tsample}_vs_{nsample}_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_maf_civic/{tsample}_vs_{nsample}.maf" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_maf_civic/{tsample}_vs_{nsample}.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_maf_civic/{tsample}_vs_{nsample}.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            code_dir=config["params"]["civic"]["code_dir"],
            category="mut",
            a_option=lambda wildcards, input: "-a %s" % input.table_alt
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_civic_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -m {params.code_dir} \
                -n {input.civic} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """
else:
    # prepare a table for all pairs tsample_vs_nsample
    rule somatic_maf_civic:
        input:
            table_alt=expand("%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
                      get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na),
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["civic"]["gene_list"],
            civic=config["params"]["civic"]["evidences"],
            rules=config["params"]["civic"]["rules_clean"]
        output:
            table_pre=temp("%s/annotation/somatic_maf_civic/all_samples_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_maf_civic/all_samples_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_maf_civic/all_samples.maf" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_maf_civic/all_samples.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_maf_civic/all_samples.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            code_dir=config["params"]["civic"]["code_dir"],
            category="mut",
            a_option=lambda wildcards, input: "-a " + " -a ".join(input.table_alt)
        threads: 1
        resources:
            queue="shortq",
            mem_mb=24000,
            time_min=90
        shell:
            """
            bash workflow/scripts/04.3_civic_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -m {params.code_dir} \
                -n {input.civic} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """


# Annnotate MAF using oncokb annotator
if config["params"]["oncokb"]["run_per_sample"]["maf"]:
    # prepare a table for each pair tsample_vs_nsample
    rule somatic_maf_oncokb:
        input:
            table_alt="%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["oncokb"]["gene_list"],
            rules=config["params"]["oncokb"]["rules_clean"]
        output:
            table_alt_pre=temp("%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}_alt_pre.tsv" % R_FOLDER),
            table_cln_pre=temp("%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}_cln_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}.maf" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_maf_oncokb/{tsample}_vs_{nsample}.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            token=config["params"]["oncokb"]["token"],
            code_dir=config["params"]["oncokb"]["code_dir"],
            category="mut",
            a_option=lambda wildcards, input: "-a %s" % input.table_alt
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_oncokb_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_alt_pre} \
                -g {output.table_cln_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -k {params.token} \
                -m {params.code_dir} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """
else:
    # prepare a table for all pairs tsample_vs_nsample
    rule somatic_maf_oncokb:
        input:
            table_alt=expand("%s/annotation/somatic_maf/{tsample}_vs_{nsample}.maf" % R_FOLDER,
                      get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na),
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["oncokb"]["gene_list"],
            rules=config["params"]["oncokb"]["rules_clean"]
        output:
            table_alt_pre=temp("%s/annotation/somatic_maf_oncokb/all_samples_alt_pre.tsv" % R_FOLDER),
            table_cln_pre=temp("%s/annotation/somatic_maf_oncokb/all_samples_cln_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_maf_oncokb/all_samples_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_maf_oncokb/all_samples.maf" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_maf_oncokb/all_samples.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_maf_oncokb/all_samples.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            token=config["params"]["oncokb"]["token"],
            code_dir=config["params"]["oncokb"]["code_dir"],
            category="mut",
            a_option=lambda wildcards, input: "-a " + " -a ".join(input.table_alt)
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_oncokb_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_alt_pre} \
                -g {output.table_cln_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -k {params.token} \
                -m {params.code_dir} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """

####
#### Copy number variants ####
####

# Annnotate SCNAs using (in-house) civic annotator
if config["params"]["civic"]["run_per_sample"]["cna"]:
    # prepare a table for each pair tsample_vs_nsample
    rule somatic_cna_civic:
        input:
            table_alt="%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["civic"]["gene_list"],
            civic=config["params"]["civic"]["evidences"],
            rules=config["params"]["civic"]["rules_clean"]
        output:
            table_pre=temp("%s/annotation/somatic_cna_civic/{tsample}_vs_{nsample}_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_cna_civic/{tsample}_vs_{nsample}_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_cna_civic/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_cna_civic/{tsample}_vs_{nsample}.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_cna_civic/{tsample}_vs_{nsample}.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            code_dir=config["params"]["civic"]["code_dir"],
            category="cna",
            a_option=lambda wildcards, input: "-a %s" % input.table_alt
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_civic_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -m {params.code_dir} \
                -n {input.civic} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """
else:
    # prepare a table for all pairs tsample_vs_nsample
    rule somatic_cna_civic:
        input:
            table_alt=expand("%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
                      get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na),
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["civic"]["gene_list"],
            civic=config["params"]["civic"]["evidences"],
            rules=config["params"]["civic"]["rules_clean"]
        output:
            table_pre=temp("%s/annotation/somatic_cna_civic/all_samples_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_cna_civic/all_samples_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_cna_civic/all_samples.tsv" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_cna_civic/all_samples.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_cna_civic/all_samples.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            code_dir=config["params"]["civic"]["code_dir"],
            category="cna",
            a_option=lambda wildcards, input: "-a " + " -a ".join(input.table_alt)
        threads: 1
        resources:
            queue="shortq",
            mem_mb=24000,
            time_min=90
        shell:
            """
            bash workflow/scripts/04.3_civic_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -m {params.code_dir} \
                -n {input.civic} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """


# Annnotate SCNAs using oncokb annotator
if config["params"]["oncokb"]["run_per_sample"]["cna"]:
    # prepare a table for each pair tsample_vs_nsample
    rule somatic_cna_oncokb:
        input:
            table_alt="%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["oncokb"]["gene_list"],
            rules=config["params"]["oncokb"]["rules_clean"]
        output:
            table_alt_pre=temp("%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}_alt_pre.tsv" % R_FOLDER),
            table_cln_pre=temp("%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}_cln_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}.tsv" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_cna_oncokb/{tsample}_vs_{nsample}.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            token=config["params"]["oncokb"]["token"],
            code_dir=config["params"]["oncokb"]["code_dir"],
            category="cna",
            a_option=lambda wildcards, input: "-a %s" % input.table_alt
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_oncokb_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_alt_pre} \
                -g {output.table_cln_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -k {params.token} \
                -m {params.code_dir} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """
else:
    # prepare a table for all pairs tsample_vs_nsample
    rule somatic_cna_oncokb:
        input:
            table_alt=expand("%s/calling/somatic_cnv_gene_calls_filtered/{tsample}_vs_{nsample}.tsv.gz" % R_FOLDER,
                      get_allowed_pairs_tumor_normal(), tsample=tsamples, nsample=nsamples_na),
            table_cln="config/tumor_normal_pairs.tsv",
            table_gen=config["params"]["oncokb"]["gene_list"],
            rules=config["params"]["oncokb"]["rules_clean"]
        output:
            table_alt_pre=temp("%s/annotation/somatic_cna_oncokb/all_samples_alt_pre.tsv" % R_FOLDER),
            table_cln_pre=temp("%s/annotation/somatic_cna_oncokb/all_samples_cln_pre.tsv" % R_FOLDER),
            table_run=temp("%s/annotation/somatic_cna_oncokb/all_samples_run.tsv" % R_FOLDER),
            table_pos="%s/annotation/somatic_cna_oncokb/all_samples.tsv" % R_FOLDER,
        benchmark:
            "%s/annotation/somatic_cna_oncokb/all_samples.tsv" % B_FOLDER
        log:
            "%s/annotation/somatic_cna_oncokb/all_samples.log" % L_FOLDER
        conda:
            "../envs/python.yaml"
        params:
            token=config["params"]["oncokb"]["token"],
            code_dir=config["params"]["oncokb"]["code_dir"],
            category="cna",
            a_option=lambda wildcards, input: "-a " + " -a ".join(input.table_alt)
        threads: 1
        resources:
            queue="shortq",
            mem_mb=4000,
            time_min=20
        shell:
            """
            bash workflow/scripts/04.3_oncokb_annotate.sh \
                {params.a_option} \
                -b {input.table_cln} \
                -c {input.table_gen} \
                -d {output.table_alt_pre} \
                -g {output.table_cln_pre} \
                -e {output.table_run} \
                -f {output.table_pos} \
                -k {params.token} \
                -m {params.code_dir} \
                -o {input.rules} \
                -t {params.category} \
                -l {log}
            """
