# @created: 20 Aug 21
# @modified: 20 Aug 21
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rprism))

source("workflow/functions/io.R")

# functions ============================================================================================================

prepare_burden_table <- function(df_fus, df_sam, df_cln, col_subject_id, col_tt, tt_keep){
  # count fusions per tumor
  df_burden <- df_fus %>%
    group_by(.data[[col_subject_id]], .data[[col_tt]]) %>%
    summarize(Fusion_Burden=n(), .groups="keep")

  # add samples with no detected fusions
  df_sam <- left_join(df_sam, df_cln[,c(col_subject_id, col_tt, "Sample_Id_RNA_T")], by=c("Sample_Id"="Sample_Id_RNA_T"))
  df_burden_mis <- df_sam %>%
    select(all_of(c(col_subject_id, col_tt))) %>%
    filter(!.data[[col_subject_id]] %in% df_burden[[col_subject_id]]) %>%
    mutate(Fusion_Burden=0)

  rbind(df_burden, df_burden_mis)
}


prepare_burden_tables <- function(dfs_fus, dfs_sam, dfs_cln, col_subject_id, col_tt, tt_keep, analysis){
  dfs_burden <- list()

  for (cohort in names(dfs_fus)){
    df_burden <- prepare_burden_table(df_fus=dfs_fus[[cohort]],
                                      df_sam=dfs_sam[[cohort]],
                                      df_cln=dfs_cln[[cohort]],
                                      col_subject_id=col_subject_id,
                                      col_tt=col_tt,
                                      tt_keep=tt_keep)

    dfs_burden[[cohort]] <- df_burden
  }

  dfs_burden
}


main <- function(args){
  col_tt <- "Tumor_Type"

  # load input counts
  df_counts_tt <- read_tsv(args$input_counts, col_types=cols())
  if (args$use_all_samples){
    tt_keep <- df_counts_tt %>% pull(.data[[col_tt]])
  } else {
    tt_keep <- df_counts_tt %>% filter(Use_For_burden==1) %>% pull(.data[[col_tt]])
  }

  # load fusions
  if (args$use_all_samples){
    dfs_sam <- load_samples_for_analysis(args$input_cohorts, args$input_samples, NULL)
    dfs_fus <- load_fusions_for_analysis(args$input_cohorts, args$input_fusions, dfs_sam, "burden")
  } else {
    dfs_sam <- load_samples_for_analysis(args$input_cohorts, args$input_samples, "burden")
    dfs_fus <- load_fusions_for_analysis(args$input_cohorts, args$input_fusions, args$input_samples, "burden")
  }

  # prepare for plot
  dfs_cln <- setNames(lapply(args$input_cohorts, load_cln), args$input_cohorts)
  dfs_cln <- lapply(dfs_cln, function(df) df %>% rename(Tumor_Type=Project_TCGA_More))
  dfs_burden <- prepare_burden_tables(dfs_fus=dfs_fus, dfs_sam=dfs_sam, dfs_cln=dfs_cln, col_subject_id="Subject_Id",
                                      col_tt=col_tt, tt_keep=tt_keep)

  # save
  for (i in seq(1,length(args$input_cohorts))){
    write.table(dfs_burden[[args$input_cohorts[i]]], args$output[i], quote=F, row.names=F, sep="\t")
    cat(paste("-file saved at", args$output[i], "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Prepare tables for the fusions burden plots.')
  parser$add_argument("--input_cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--input_samples", nargs="+", help="Paths to input sample tables.",
                      default=c("workflow/results/tables/prism/selection_samples.tsv",
                                "workflow/results/tables/met500/selection_samples.tsv",
                                "workflow/results/tables/tcga/selection_samples.tsv"))
  parser$add_argument("--input_fusions", nargs="+", help="Paths to input fusion tables.",
                      default=c("workflow/results/tables/prism/fusions_preprocessed.tsv",
                                "workflow/results/tables/met500/fusions_preprocessed.tsv",
                                "workflow/results/tables/tcga/fusions_preprocessed.tsv"))
  parser$add_argument("--input_counts", type="character", help="Path to count table of tumor types.",
                      default="workflow/results/tables/selection/selection_tumor_types.tsv")
  parser$add_argument("--use_all_samples", action="store_true", default=F, help="If specified, all samples are used.")
  parser$add_argument("--output", nargs="+", help="Path to where output tables will be saved.",
                      default=c("workflow/results/tables/prism/burden_table_sel.tsv",
                                "workflow/results/tables/met500/burden_table_sel.tsv",
                                "workflow/results/tables/tcga/burden_table_sel.tsv"))
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
