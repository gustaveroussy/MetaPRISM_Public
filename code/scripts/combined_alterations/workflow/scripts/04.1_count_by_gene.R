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

  # for alterations that are not already equal to gene symbol, add gene symbol as prefix
  mask <- df_alt[[col_alt]] != df_alt[[col_gen]]
  df_alt[mask,col_alt] <- df_alt %>% filter(mask) %>%
    unite("Unite", all_of(c(col_gen, col_alt)), sep=" ") %>% pull(var="Unite")

  # simplifcation rules
  # a. if for a gene XXXX, `col_alt` is only "XXXX Oncogenic, no Level", the latter is simplified to "XXXX"
  # b. if for a gene XXXX, `col_alt` is either "XXXX" or "XXXX Oncogenic, no Level", the latter is simplified to "XXXX"
  for (gen in unique(df_alt[[col_gen]])){
    mask_gen <- df_alt[[col_gen]]==gen
    df_alt_gen <- df_alt %>% filter(mask_gen) %>% select(all_of(c(col_gen, col_alt))) %>% distinct()
    alt_gen <- df_alt_gen %>% pull(var=col_alt)
    gen_no_lvl <- paste(gen, "Oncogenic, no Level")
    if (all(alt_gen %in% c(gen, gen_no_lvl))){
      df_alt <- df_alt %>% mutate(!!col_alt:=ifelse(mask_gen, gen, .data[[col_alt]]))
    }
  }

  # get counts per each value of "col_alt"
  df_cnt <- df_alt %>% group_by(Subject_Id, .data[[col_alt]]) %>% summarize(Count=n(), .groups="drop")  %>%
    spread(.data[[col_alt]], Count)

  # add subjects with no events
  df_cnt <- bind_rows(df_cnt, tibble(Subject_Id=setdiff(df_sam$Subject_Id, df_cnt$Subject_Id)))

  # replace NA by 0
  df_cnt <- df_cnt %>% replace(is.na(.), 0)

  # save
  write.table(df_cnt, args$output, sep="\t", quote=F, row.names=F)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build table of alteration counts per gene.')
  parser$add_argument("--alterations", type="character", help="Path to mutations table.",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations.tsv")
  parser$add_argument("--samples", type="character", help="Path to table of selected samples.",
                      default="../../../results/combined_alterations/selection/selection_samples_prism.tsv")
  parser$add_argument("--data_type", type="character", help="Type of data to be used. Choose 'DNA' or 'DNA_RNA'",
                      default="DNA")
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--output", type="character", help="Path to output left bar data file.",
                      default="../../../results/combined_alterations/count/count_by_gene_DNA_annotated_prism.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
