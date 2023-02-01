# @created: 19 Apr 22
# @modified: 27 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

add_counts_total <- function(df_counts, df_rna){
  df_counts_total <- df_rna %>% select(Sample_Id, Fusion_Id) %>% distinct() %>% group_by(Sample_Id) %>%
    summarize(Count=n())
  df_counts_total <- df_counts_total %>% rename(Total:=Count)
  df_counts <- left_join(df_counts, df_counts_total, how="left", 
                         by=c("Sample_Id"))

  df_counts
}

main <- function(args){
  df_sam <- load_table(args$samples)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
    cat(paste("-number of samples that pass QC:", nrow(df_sam), "\n"))
  }
  df_rna <- load_table(args$fusions)

  # select samples in samples table
  df_rna <- df_rna %>% filter(Sample_Id %in% df_sam$Sample_Id)
  cat(paste("-number of variants:", nrow(df_rna), "\n"))

  # select samples in cln table
  df_cln <- load_cln(study=args$cohort)
  df_rna <- df_rna %>% filter(Sample_Id %in% df_cln$Sample_Id_RNA_T)
  df_sam <- df_sam %>% filter(Sample_Id %in% df_cln$Sample_Id_RNA_T)

  # get counts table
  df_counts <- df_sam %>% select(Sample_Id)
  df_counts <- add_counts_total(df_counts, df_rna)
  df_counts <- df_counts %>% replace(is.na(.), 0)

  # remove duplicates if any (TCGA)
  df_counts <- df_counts %>% group_by(Sample_Id) %>% summarize(Total=max(Total))

  # save
  write.table(df_counts, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {

  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="tcga")
  parser$add_argument("--samples", type="character", help="Path to sample list table.",
                      default="../../../data/tcga/rna/fusions/sample_list.tsv")
  parser$add_argument("--fusions", type="character", help="Path to fusions table.",
                      default="../../../data/tcga/rna/fusions/tcga_annotated_filtered.tsv.gz")
  parser$add_argument("--output", type="character", help="Path to output counts table.",
                      default="../../../results/fusions_analysis/count/count_total_RNA_all_tcga.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
