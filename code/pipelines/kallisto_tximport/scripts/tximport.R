#!/bin/R

# Loading libraries
# See https://github.com/gevaertlab/RNASeq_pipeline/blob/master/scripts/tximport.R#2
base::library(package="tximport", quietly=TRUE);


# Loading tx2gene table
# See https://github.com/gevaertlab/RNASeq_pipeline/blob/master/scripts/tximport.R#20
tx2gene <- utils::read.table(
  file = base::as.character(x = snakemake@input[["tx2gene"]]),
  header = FALSE,
  sep = "\t"
);


# Creating tximport object
# See https://github.com/gevaertlab/RNASeq_pipeline/blob/master/scripts/tximport.R#21
kallisto_files <- sapply(
    snakemake@input[["quant"]],
    function(path) base::as.character(x = path)
);

kallisto <- tximport::tximport(
  files = kallisto_files,
  type = "kallisto",
  tx2gene = tx2gene,
  ignoreAfterBar = TRUE
);


# Save results
# See https://github.com/gevaertlab/RNASeq_pipeline/blob/master/scripts/tximport.R#24
utils::write.csv(
  kallisto,
  file=snakemake@output[["txt"]],
  row.names = T
);

base::saveRDS(
  kallisto,
  file=snakemake@output[["txi"]],
)
