rule logistic_regression:
    benchmark:
        "%s/logistic_regression.csv" % B_FOLDER
    log:
        "%s/logistic_regression.log" % L_FOLDER
    input:
        cln = "%s/prism/clinical/curated/cln_prism_in_design_curated.tsv" % D_FOLDER,
        sig = "%s/projection_known_signatures/MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_prism.tsv" % R_FOLDER,
        env = "../common/logs/setup_conda.done"
    conda:
        config["setup"]["MetaPrism"]
    output:
        out = "%s/associations/logistic_regression_platinum.pdf" % R_FOLDER,
        out_paper = "%s/F2b_bot.pdf" % F_FOLDER
    threads: 1
    shell:
        """Rscript workflow/scripts/05_logistic_regression.R \
            --cln {input.cln} \
            --sig {input.sig} \
            --output {output.out} \
            --output_paper {output.out_paper} \
            --log {log}"""

