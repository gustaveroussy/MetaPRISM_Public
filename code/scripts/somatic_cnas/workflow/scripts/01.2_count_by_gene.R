# @modified: 06 Oct 21
# @modified: 06 Oct 21
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================


get_corrected_sample_id <- function(sample_id){
  sample_id <- gsub("^Sample_","",sample_id)
  sample_id <- gsub("_run2$","",sample_id)

  if (sample_id=="M727-T1-ARN"){
    return("M727RE-T1-ARN")
  } else if (sample_id=="MR012-T-ARN"){
    return("MR012-T2-ARN")
  } else if (sample_id=="MR254_T2-ADN"){
    return("MR254-T2-ADN")
  } else {
    return(sample_id)
  }
}


add_counts_by_gene <- function(df_counts, df_cna, df_genes){
  gene_names <- sort(unique(df_genes$Hugo_Symbol))
  df_counts_gene <- df_cna %>% filter(HUGO_SYMBOL %in% gene_names) %>% distinct(SAMPLE_ID, HUGO_SYMBOL, ALTERATION) %>% 
    mutate(Count=ifelse(ALTERATION=="Deletion", -1, 1)) %>% select(-ALTERATION) %>% spread(HUGO_SYMBOL, Count)
  df_counts <- left_join(df_counts, df_counts_gene, how="left", by=c("Sample_Id"="SAMPLE_ID"))

  df_counts
}


add_counts_total <- function(df_counts, df_cna){
  df_counts_total <- df_cna %>% group_by(SAMPLE_ID, ALTERATION) %>% summarize(Count=n()) %>% spread(ALTERATION, Count)
  df_counts_total <- df_counts_total %>% rename(Total_Amplication=Amplification, Total_Deletion=Deletion)
  df_counts <- left_join(df_counts, df_counts_total, how="left", by=c("Sample_Id"="SAMPLE_ID"))

  df_counts
}


main <- function(args){
  df_cna <- load_table(args$cnas)
  df_cna$SAMPLE_ID <- sapply(df_cna$SAMPLE_ID, get_corrected_sample_id, USE.NAMES=F)
  df_genes <- load_table(args$genes) %>% filter(Inclusion_Level==1)
  df_samples <- load_table(args$samples)

  df_counts <- df_samples[,c("Sample_Id")]
  df_counts <- add_counts_by_gene(df_counts, df_cna, df_genes)
  df_counts <- add_counts_total(df_counts, df_cna)
  df_counts <- df_counts %>% replace(is.na(.), 0)

  # save
  write.table(df_counts, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {

  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cnas", type="character", help="Path to cnas table.",
                      default="../../../data/prism/wes/facets/prism_facets_genes_oncokb_civic.xlsx")
  parser$add_argument("--samples", type="character", help="Path to table of selected samples.",
                      default="../../../results/somatic_cnas/selection/selection_samples_prism.tsv")
  parser$add_argument("--genes", type="character", help="Path to input cnas table.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_cnas/count/count_by_gene_oncokb_civic_oncokb_cosmic.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
