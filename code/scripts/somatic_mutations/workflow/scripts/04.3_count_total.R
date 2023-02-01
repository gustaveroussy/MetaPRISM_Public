# @created: 03 Nov 21
# @modified: 26 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

add_counts_total <- function(df_counts, df_maf){
  df_counts_total <- df_maf %>% group_by(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% summarize(Count=n())
  df_counts_total <- df_counts_total %>% rename(Total:=Count)
  df_counts <- left_join(df_counts, df_counts_total, how="left", 
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
  df_counts <- df_sam %>% select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode)
  df_counts <- add_counts_total(df_counts, df_maf)
  df_counts <- df_counts %>% replace(is.na(.), 0)

  # save
  write.table(df_counts, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {

  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="tcga")
  parser$add_argument("--samples", type="character", help="Path to sample list table.",
                      default="../../../data/tcga/wes/somatic_maf/sample_list.tsv")
  parser$add_argument("--mutations", type="character", help="Path to mutations table.",
                      default="../../../data/tcga/wes/somatic_maf/somatic_calls.maf.gz")
  parser$add_argument("--selection", type="character", default="all", help="For selecting mutations.")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_mutations/count/count_total_DNA_all_tcga.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
