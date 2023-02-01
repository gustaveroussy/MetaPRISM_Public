rule setup_conda:
    conda:
        "../envs/MetaPrism.yaml"
    params:
        rprism = "../../functions/rprism/rprism_2.2.5.tar.gz",
        rprismtools = "../../functions/rprismtools/rprismtools_1.5.0.tar.gz",
        forester = "../../functions/tools/forester/forester_0.3.0.tar.gz",
        tableextra = "../../functions/tools/tableExtra",
        pyprism = "../../functions/pyprism",
        pyprismtools = "../../functions/pyprismtools",
        prettypy = "../../functions/tools/PrettyPy",
        comut = "../../functions/tools/comut",
    output:
        touch("../common/logs/setup_conda.done")
    shell:
        """
        R CMD INSTALL {params.rprism}
        R CMD INSTALL {params.rprismtools}
        R -e 'install.packages("robvis", repo="https://cran.biotools.fr")'
        R CMD INSTALL {params.forester}
        R CMD INSTALL {params.tableextra}
        Rscript -e 'install.packages("svglite", repos="https://cran.biotools.fr")'
        pip install -e {params.pyprism}
        pip install -e {params.pyprismtools}
        pip install -e {params.prettypy}
        pip install -e {params.comut}
        """


rule setup_conda_signatures:
    conda:
        "../envs/Signatures.yaml"
    params:
        rprism = "../../functions/rprism/rprism_2.2.5.tar.gz",
        ultisig = "../../functions/tools/UltiSig/UltiSig_1.0.3.tar.gz",
        sigprofilerjulia = "../mutational_signatures/external/sigprofilerjulia"
    output:
        touch("../common/logs/setup_conda_signatures.done")
    shell:
        """
        R CMD INSTALL {params.rprism}
        Rscript -e 'install.packages("coneproj", repos="https://cran.biotools.fr")'
        R CMD INSTALL {params.ultisig}
        julia {params.sigprofilerjulia}/julia_dependencies.jl
        """
