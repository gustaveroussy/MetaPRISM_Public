def get_threads_projection(w):
    if w.pkg=="deconstructSigs":
        return 20
    else:
        return 1


rule run_projection:
    wildcard_constraints:
        basis = "|".join([re.escape(x) for x in config["projection"]["bases"]]),
        cmode = "|".join([re.escape(x) for x in config["projection"]["cmodes"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]]),
        mmode = "|".join([re.escape(x) for x in config["projection"]["mmodes"]])
    benchmark:
        "%s/projection_{pkg}_{basis}_{mmode}_{cmode}_{cohort}.csv" % B_FOLDER
    log:
        "%s/projection_{pkg}_{basis}_{mmode}_{cmode}_{cohort}.log" % L_FOLDER
    input:
        count_mut = "%s/counts_mutations/counts_mutations_{mmode}_{cmode}_{cohort}.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda_signatures.done"
    conda:
        config["setup"]["Signatures"]
    output:
        R = "%s/projection_known_signatures/{pkg}/counts_mutations_{basis}_{mmode}_{cmode}_{cohort}.tsv" % R_FOLDER,
        H = "%s/projection_known_signatures/{pkg}/counts_signatures_{basis}_{mmode}_{cmode}_{cohort}.tsv" % R_FOLDER
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "02:00:00"
    threads: get_threads_projection
    shell:
        """Rscript workflow/scripts/02_run_projection.R \
            --count {input.count_mut} \
            -b {wildcards.basis} \
            -k {wildcards.pkg} \
            --n_cores {threads} \
            --filepath_H {output.H} \
            --filepath_R {output.R} \
            --log {log}"""


rule sparsify_w_sigprofilerjulia:
    wildcard_constraints:
        pkg = "|".join([re.escape(x) for x in config["projection"]["pkgs"]]),
        basis = "|".join([re.escape(x) for x in config["projection"]["bases"]]),
        cmode = "|".join([re.escape(x) for x in config["projection"]["cmodes"]]),
        cohort = "|".join([re.escape(x) for x in config["data"]["cohorts"]]),
        mmode = "|".join([re.escape(x) for x in config["projection"]["mmodes"]])
    benchmark:
        "%s/sparsify_w_sigprofilerjulia_{pkg}_cosmic_sbs_96_v3.2_sbs_96_min_mut_{cohort}.csv" % B_FOLDER
    log:
        "%s/sparsify_w_sigprofilerjulia_{pkg}_cosmic_sbs_96_v3.2_sbs_96_min_mut_{cohort}.log" % L_FOLDER
    input:
        count_mut = "%s/counts_mutations/counts_mutations_sbs_96_min_mut_{cohort}.tsv" % R_FOLDER,
        count_sig = "%s/projection_known_signatures/{pkg}/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_{cohort}.tsv" % R_FOLDER,
        basis = config["data"]["resources"]["bases"]["cosmic_sbs_96_v3.2"],
        env = "../common/logs/setup_conda_signatures.done"
    conda:
        config["setup"]["Signatures"]
    output:
        "%s/projection_known_signatures/{pkg}/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_{cohort}.tsv" % R_FOLDER
    params:
        sigprofilerjulia = config["projection"]["sigprofilerjulia"]
    resources:
        partition = "cpu_med",
        mem_mb = 24000,
        time = "03:00:00"
    threads: 20
    shell:
        """julia --procs={threads} workflow/scripts/03_sparsify_w_sigprofilerjulia.jl \
            --count_mut {input.count_mut} \
            --count_sig {input.count_sig} \
            --basis {input.basis} \
            --sigprofilerjulia {params.sigprofilerjulia} \
            --output {output} &> {log}"""


rule table_extra:
    benchmark:
        "%s/table_extra_MutationalPatterns_cosmic_sbs_96_v3.2_sbs_96_min_mut.csv" % B_FOLDER
    log:
        "%s/table_extra_MutationalPatterns_cosmic_sbs_96_v3.2_sbs_96_min_mut.log" % L_FOLDER
    input:
        clns = expand("%s/{cohort}/clinical/curated/cln_{cohort}_in_design_curated.tsv" % D_FOLDER,
            cohort=["prism", "met500", "tcga"]),
        sigs = expand("%s/projection_known_signatures/MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_{cohort}.tsv" % R_FOLDER,
            cohort=["prism", "met500", "tcga"]),
        sams = expand("%s/selection/selection_samples_{cohort}.tsv" % R_FOLDER,
            cohort=["prism", "met500", "tcga"]),
        counts = "%s/selection/selection_tumor_types.tsv" % R_FOLDER,
        target_bed = "%s/resources/target_files/all_targets_intersect_padded_10n.bed" % D_FOLDER,
        aetiologies = config["data"]["resources"]["aetiologies"]["cosmic_sbs_96_v3.2"],
    conda:
        config["setup"]["MetaPrism"]
    output:
        outs = expand("%s/projection_known_signatures/MutationalPatterns/table_extra_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_{table}.pdf" % R_FOLDER, table=["full","plat"]),
        outs_paper = expand("%s/F{x}.pdf" % F_FOLDER, x=["S4a", "S4b", "2b_top"])
    params:
        cohorts = ["prism", "met500", "tcga"],
    threads: 1
    shell:
        """Rscript workflow/scripts/04_table_extra.R \
            --clns {input.clns} \
            --sigs {input.sigs} \
            --sams {input.sams} \
            --counts {input.counts} \
            --target_bed {input.target_bed} \
            --aetiologies {input.aetiologies} \
            --outputs {output.outs} \
            --outputs_paper {output.outs_paper} \
            --log {log}"""
