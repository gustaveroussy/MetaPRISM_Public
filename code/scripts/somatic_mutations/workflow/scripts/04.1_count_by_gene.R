# @created: 07 Oct 21
# @modified: 26 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

add_counts_by_gene <- function(df_counts, df_maf, df_genes){
  gene_names <- sort(unique(df_genes$Hugo_Symbol))
  df_counts_gene <- df_maf %>% filter(Hugo_Symbol %in% gene_names) %>% 
    distinct(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Hugo_Symbol) %>% 
    mutate(Count=1) %>% spread(Hugo_Symbol, Count)
  df_counts <- left_join(df_counts, df_counts_gene, how="left",
                         by=c("Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode"))

  df_counts
}

main <- function(args){
  df_sam <- load_table(args$samples)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
    cat(paste("-number of samples that pass QC:", nrow(df_sam), "\n"))
  }
  df_sam <- df_sam %>% rename(Tumor_Sample_Barcode=Tumor_Sample_Id, Matched_Norm_Sample_Barcode=Normal_Sample_Id)
  df_maf <- load_table(args$mutations, header_prefix="##")

  # select samples in samples table
  df_maf <- df_maf %>% unite("DNA_P", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)
  df_sam <- df_sam %>% unite("DNA_P", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)
  df_maf <- df_maf %>% filter(DNA_P %in% df_sam$DNA_P) %>% select(-DNA_P)
  cat(paste("-number of variants:", nrow(df_maf), "\n"))

  # select samples in clinical table
  df_sam <- preprocess_wes_mut(df_sam, args$cohort, select_pairs=T, verbose=T)
  df_maf <- preprocess_wes_mut(df_maf, args$cohort, select_pairs=T, selection_mut=args$selection, verbose=T)

  # get counts table
  df_genes <- load_table(args$genes) %>% filter(Inclusion_Level==1)
  df_counts <- df_sam %>% select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode)
  df_counts <- add_counts_by_gene(df_counts, df_maf, df_genes)
  df_counts <- df_counts %>% replace(is.na(.), 0)

  # save
  write.table(df_counts, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build table of mutation counts per gene.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--samples", type="character", help="Path to sample list table.",
                      default="../../../data/prism/wes/somatic_maf/sample_list.tsv")
  parser$add_argument("--mutations", type="character", help="Path to mutations table.",
                      default="../../../data/prism/wes/somatic_maf/somatic_calls.maf.gz")
  parser$add_argument("--selection", type="character", default="non_synonymous", help="For selecting mutations.")
  parser$add_argument("--genes", type="character", help="Path to genes table.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_mutations/count/count_by_gene_oncokb_cosmic_DNA_all_prism.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
