# @created: 07 Oct 21
# @modified: 06 Jul 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

select_samples <- function(df_agg, df_sam, col_use){
  cat("-selecting samples ...\n")
  df_sam <- df_sam %>% filter(.data[[col_use]]==1)
  n_row_bef <- nrow(df_agg)
  df_agg <- df_agg %>% filter(Sample_Id %in% df_sam$Sample_Id)
  n_row_aft <- nrow(df_agg)
  cat("-INFO: selected", paste0(n_row_aft, "/", n_row_bef), "lines by selecting samples\n")
  cat("-INFO: number of unique subjects after selection:", df_agg %>% distinct(Subject_Id) %>% nrow(), "\n")
  df_agg
}


main <- function(args){
  # load tables
  df_alt <- load_table(args$alterations)
  df_sam <- load_table(args$samples)

  # select cohort
  df_alt <- df_alt %>% filter(Cohort==args$cohort)

  # select data type and samples
  if (args$data_type=="DNA"){
    df_alt <- df_alt %>% filter(Alteration_Category %in% c("Amplification", "Deletion", "Ins", "Del", "Mut"))
    df_alt <- select_samples(df_alt, df_sam, col_use="Use_heatmap_dna")
    df_sam <- df_sam %>% filter(Use_heatmap_dna==1)
  } else if (args$data_type=="DNA_RNA"){
    df_alt <- df_alt %>% filter(!Alteration_Category %in% c("MSI", "TMB"))
    df_alt <- select_samples(df_alt, df_sam, col_use="Use_heatmap_all")
    df_sam <- df_sam %>% filter(Use_heatmap_all==1)
  } else {
    stop("-ERROR: only DNA and DNA_RNA are supported for the argument --data_type")
  }

  # define util variable names
  col_lvl <- "Sen_Level_Simple"
  col_gen <- paste0("Hugo_Symbol_", col_lvl)
  col_alt <- paste0("Alteration_", col_lvl)
  col_tt <- "Tumor_Type"
  col_id <- "Subject_Id"

  # add columns discriminating between Annotated_Level and Annotated_Tier1/Tier2/Tier3
  df_alt <- df_alt %>% mutate(Level=ifelse(is.na(.data[[col_lvl]]), "Annotated_No_Level",
                                           paste0("Annotated_", .data[[col_lvl]])))

  # get counts for each pathway
  df_pathways <- load_table(args$pathways)
  pathway_names <- sort(unique(df_pathways$Pathway_Name))

  df_cnt <- df_sam[,c("Subject_Id")] %>% distinct()
  for (pathway_name in pathway_names){
    members <- df_pathways %>% filter(Pathway_Name==pathway_name) %>% pull(Symbol)
    df_cnt_pathway <- df_alt %>% filter(.data[[col_gen]] %in% members) %>% group_by(Subject_Id, Level) %>%
      summarize(Count=n(), .groups="drop") %>% spread(Level, Count)
    colnames(df_cnt_pathway) <- c("Subject_Id", 
                                  paste0(pathway_name, "_", colnames(df_cnt_pathway)[2:ncol(df_cnt_pathway)]))
    df_cnt <- left_join(df_cnt, df_cnt_pathway, by=c("Subject_Id"))
  }

  # replace NA by 0
  df_cnt <- df_cnt %>% replace(is.na(.), 0)

  # save
  write.table(df_cnt, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build table of alterations counts per pathway.')
  parser$add_argument("--alterations", type="character", help="Path to mutations table.",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations.tsv")
  parser$add_argument("--samples", type="character", help="Path to table of selected samples.",
                      default="../../../results/combined_alterations/selection/selection_samples_tcga.tsv")
  parser$add_argument("--data_type", type="character", help="Type of data to be used. Choose 'DNA' or 'DNA_RNA'",
                      default="DNA_RNA")
  parser$add_argument("--pathways", type="character", help="Path to pathways table.",
                      default="../../../data/resources/pathways/data/sanchez-vega_cell_2018_curated.tsv")
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="tcga")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
            default="../../../results/combined_alterations/count/count_by_pathway_sanchez_vega_DNA_annotated_tcga.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}

