# @created: 03 Nov 21
# @modified: 03 Nov 21
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))


# functions ============================================================================================================

add_counts_total <- function(df_counts, df_maf){
  df_counts_total <- df_maf %>% group_by(Tumor_Sample_Barcode) %>% summarize(Count=n())
  df_counts_total <- df_counts_total %>% rename(Total:=Count)
  df_counts <- left_join(df_counts, df_counts_total, how="left", by=c("Sample_Id"="Tumor_Sample_Barcode"))

  df_counts
}

main <- function(args){
  df_maf <- load_table(args$mutations, header_prefix="##", guess_max=1e5)
  df_maf <- preprocess_wes_mut(df_maf, args$cohort, select_pairs=T, selection_mut=args$selection_mut, verbose=T)
  df_samples <- load_table(args$samples)

  df_counts <- df_samples[,c("Sample_Id")]
  df_counts <- add_counts_total(df_counts, df_maf)
  df_counts <- df_counts %>% replace(is.na(.), 0)

  # save
  write.table(df_counts, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {

  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--mutations", type="character", help="Path to mutations table.",
                      default="../../../data/tcga/wes/oncotator/tcga_oncotator_oncokb_civic.xlsx")
  parser$add_argument("--selection_mut", type="character", default="all", help="For selecting mutations.")
  parser$add_argument("--samples", type="character", help="Path to table of selected samples.",
                      default="../../../results/somatic_mutations/selection/selection_samples_tcga.tsv")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_mutations/count/count_total_for_all_tcga.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
