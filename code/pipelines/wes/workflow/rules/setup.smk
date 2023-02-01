rule setup_r:
    conda:
        "../envs/r.yaml"
    output:
        touch("%s/setup_r.done" % L_FOLDER)
    resources:
        queue="shortq",
        mem_mb=16000,
        time_min=60
    shell:
        """
        Rscript -e 'devtools::install_github("mskcc/facets-suite")'
        """

