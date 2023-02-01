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



add_counts_by_pathway <- function(df_counts, df_cna, df_pathways){
  pathway_names <- sort(unique(df_pathways$Pathway_Name))
  for (pathway_name in pathway_names){
    members <- df_pathways %>% filter(Pathway_Name==pathway_name) %>% pull(Symbol)
    df_counts_pathway <- df_cna %>% filter(HUGO_SYMBOL %in% members) %>% group_by(SAMPLE_ID) %>%
      summarize(!!pathway_name:=n())
    df_counts <- left_join(df_counts, df_counts_pathway, how="left", by=c("Sample_Id"="SAMPLE_ID"))
  }

  df_counts
}


add_counts_total <- function(df_counts, df_cna){
  df_counts_total <- df_cna %>% group_by(SAMPLE_ID) %>% summarize(Count=n())
  df_counts_total <- df_counts_total %>% rename(Total=Count)
  df_counts <- left_join(df_counts, df_counts_total, how="left", by=c("Sample_Id"="SAMPLE_ID"))

  df_counts
}


main <- function(args){
  df_cna <- load_table(args$cnas)
  df_cna$SAMPLE_ID <- sapply(df_cna$SAMPLE_ID, get_corrected_sample_id, USE.NAMES=F)
  df_pathways <- load_table(args$pathways)
  df_samples <- load_table(args$samples)

  df_counts <- df_samples[,c("Sample_Id")]
  df_counts <- add_counts_by_pathway(df_counts, df_cna, df_pathways)
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
  parser$add_argument("--pathways", type="character", help="Path to input cnas table.",
                      default="../../../data/resources/pathways/data/sanchez-vega_cell_2018_curated.tsv")
  parser$add_argument("--prefix", type="character", help="Prefix to count columns in output file",
                      default="Somatic_Mutations_filtered")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_cnas/count/count_by_pathway_filtered_sanchez_vega_prism.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
