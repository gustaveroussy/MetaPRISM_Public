rule setup_conda:
    conda:
        config["setup"]["MetaPrism"]
    params:
        rprism = config["setup"]["local_packages"]["rprism"],
        rprismtools = config["setup"]["local_packages"]["rprismtools"],
        tableextra = config["setup"]["local_packages"]["tableextra"]
    output:
        touch("../common/logs/setup_conda.done")
    shell:
        """
        R -e 'install.packages(c("{params.rprism}", "{params.rprismtools}", "{params.tableextra}"), \
            repos=NULL, type="source")'
        R -e 'devtools::install_github("rdboyes/forester")'
        """

rule setup_conda:
    conda:
        config["setup"]["MetaPrism"]
    params:
        pyprism = config["setup"]["local_packages"]["pyprism"],
        pyprismtools = config["setup"]["local_packages"]["pyprismtools"],
        prettypy = config["setup"]["local_packages"]["prettypy"]
    output:
        touch("../common/logs/setup_conda.done")
    shell:
        """
        pip install -e {params.pyprism}
        pip install -e {params.pyprismtools}
        pip install -e {params.prettypy}
        """
