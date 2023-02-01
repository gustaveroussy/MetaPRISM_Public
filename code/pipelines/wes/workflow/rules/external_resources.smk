# Download databases required for annotating with annovar
# See https://annovar.openbioinformatics.org/en/latest/user-guide/filter
rule annovar_downdb:
    output:
        "%s/humandb/hg19_{resource}.txt" % config["params"]["annovar"]["path"],
        "%s/humandb/hg19_{resource}.txt.idx" % config["params"]["annovar"]["path"]
    benchmark:
        "%s/annotation/annovar_downdb/{resource}.tsv" % B_FOLDER
    log:
        "%s/annotation/annovar_downdb/{resource}.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        annovar_dir=config["params"]["annovar"]["path"],
        buildver=config["params"]["annovar"]["buildver"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=8000,
        time_min=90
    shell:
        """{params.annovar_dir}/annotate_variation.pl -buildver {params.buildver} \
            -downdb \
            -webfrom annovar {wildcards.resource} {params.annovar_dir}/humandb 2> {log}"""


# Install MANTIS
# See https://github.com/OSU-SRLab/MANTIS
rule install_mantis:
    input:
        ref=config["ref"]["fasta"],
        bed=config["params"]["mantis"]["ms_bed"]
    output:
        "%s/mantis.py" % config["params"]["mantis"]["code_dir"],
        "%s/%s_microsatellites.bed" % (config["params"]["mantis"]["ms_dir"], config["params"]["mantis"]["buildver"]),
        expand("%s/%s_microsatellites_{target}.bed" %
            (config["params"]["mantis"]["ms_dir"],
             config["params"]["mantis"]["buildver"]), target=get_target_names()),
        "%s/%s_microsatellites_%s" %
            (config["params"]["mantis"]["ms_dir"],
             config["params"]["mantis"]["buildver"],
             os.path.basename(config["params"]["mantis"]["ms_bed"]))
    benchmark:
        "%s/calling/install_mantis/mantis.tsv" % B_FOLDER
    log:
        "%s/calling/install_mantis/mantis.log" % L_FOLDER
    conda:
        "../envs/main.yaml"
    params:
        code_dir=config["params"]["mantis"]["code_dir"],
        ms_dir=config["params"]["mantis"]["ms_dir"],
        buildver=config["params"]["mantis"]["buildver"],
        targets_name=" -n " + " -n ".join(get_target_names()),
        targets_bed=" -b " + " -b ".join(get_target_files(get_target_names(), file="bed")),
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=30
    shell:
        """
        bash workflow/scripts/02.1_install_mantis.sh \
            -f {input.ref} \
            -c {params.code_dir} \
            -s {params.ms_dir} \
            -v {params.buildver} \
            -d {input.bed} \
            {params.targets_name} {params.targets_bed} &> {log}
        """


# Download oncokb code and required databases if not done yet.
rule install_oncokb:
    output:
        "%s/CnaAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/FusionAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/MafAnnotator.py" % config["params"]["oncokb"]["code_dir"],
        "%s/cancerGeneList_oncokb_annotated.tsv" % config["params"]["oncokb"]["data_dir"]
    benchmark:
        "%s/annotation/install_oncokb/oncokb.tsv" % B_FOLDER
    log:
        "%s/annotation/install_oncokb/oncokb.log" % L_FOLDER
    params:
        code_dir=config["params"]["oncokb"]["code_dir"],
        data_dir=config["params"]["oncokb"]["data_dir"],
        gene_list=config["params"]["oncokb"]["gene_list"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=4000,
        time_min=15
    shell:
        """
        bash workflow/scripts/02.2_install_oncokb.sh \
            -c {params.code_dir} \
            -d {params.data_dir} \
            -g {params.gene_list} &> {log}
        """


# Install VEP code
# See https://m.ensembl.org/info/docs/tools/vep/script/vep_download.html
# See http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html
# See https://githubmemory.com/repo/Ensembl/ensembl-vep/issues/930
rule install_vep:
    output:
        "%s/INSTALL.pl" % config["params"]["vep"]["path"],
        # "%s/%s" % (config["params"]["vep"]["cache"], config["params"]["vep"]["fasta"])
    benchmark:
        "%s/annotation/install_vep/vep.tsv" % B_FOLDER
    log:
        "%s/annotation/install_vep/vep.log" % L_FOLDER
    conda:
        "../envs/perl.yaml"
    params:
        assembly=config["ref"]["build"],
        species=config["ref"]["species"],
        plugins=",".join(config["params"]["vep"]["plugins"]),
        cache=config["params"]["vep"]["cache"],
        release=config["params"]["vep"]["release"],
        fasta="%s/%s.gz" % (config["params"]["vep"]["cache"], config["params"]["vep"]["fasta"]),
        path=config["params"]["vep"]["path"]
    threads: 1
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=240
    shell:
        """
        bash workflow/scripts/02.3_install_vep.sh \
            -a {params.assembly} \
            -s {params.species} \
            -g {params.plugins} \
            -c {params.cache} \
            -r {params.release} \
            -f {params.fasta} \
            -p {params.path} &> {log}
        """

# Download plugins databases required for annotating with vep plugins
# See https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
rule install_vep_plugins:
    input:
        "%s/INSTALL.pl" % config["params"]["vep"]["path"]
    output:
        config["params"]["vep"]["plugins_data"]["CADD"],
        config["params"]["vep"]["plugins_data"]["Condel"],
        config["params"]["vep"]["plugins_data"]["dbNSFP"],
        config["params"]["vep"]["plugins_data"]["LoFtool"],
    benchmark:
        "%s/annotation/install_vep_plugins/vep.tsv" % B_FOLDER
    log:
        "%s/annotation/install_vep_plugins/vep.log" % L_FOLDER
    params:
        cache=config["params"]["vep"]["cache"],
        assembly=config["ref"]["build"],
        cadd_version=config["params"]["vep"]["plugins_version"]["CADD"],
        dbnsfp_version=config["params"]["vep"]["plugins_version"]["dbNSFP"]
    threads: 6
    resources:
        queue="longq",
        mem_mb=20000,
        time_min=1440
    shell:
        """
        bash workflow/scripts/02.4_install_vep_plugins.sh \
            -c {params.cache} \
            -a {params.assembly} \
            -u {params.cadd_version} \
            -v {params.dbnsfp_version} &> {log}
        """
