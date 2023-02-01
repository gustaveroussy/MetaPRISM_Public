# @created: 29 Aug 22
# @modified: 31 Aug 22
# @authors: Yoann Pradat
#
# Upset plots showing how fusions were called.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

build_indices <- function(df, field, unique_values, verbose=F, special_values=c()){
  df_indices <- data.frame()
  for (value in unique_values){
    if (verbose){
      cat(paste("building indices for", value, "...\n"))
    }

    df_sub <- df[,field,drop=F]
    if (value %in% special_values){
      rows <- 1:nrow(df_sub)
      rows <- rows[df_sub[[field]]==value]
    } else {
      rows <- which(
        rowSums(
          `dim<-`(grepl(value, as.matrix(df_sub), fixed=TRUE), dim(df_sub))
        ) > 0
      )
    }
    df_ind <- data.frame(Index=rows)
    df_ind[,"Set"] <- value
    df_indices <- rbind(df_indices, df_ind)
  }
  df_indices
}


draw_upset_plot_fusion <- function(df_fus, output, field="Algo_Wo_Breakpoint"){
  unique_values <- sort(unique(scan(text=df_fus[[field]], what="", sep=";")))
  df_fus_indices <- build_indices(df=df_fus, field=field, unique_values=unique_values)
  m_upset <- make_upset_m(df_fus_indices, field_set="Set", field_identifier="Index")
  pdf(file=output, width=12, height=6)
  draw_upset_plot(m_upset, width_set_size=3, height_top_annot=3)
  dev.off()
  cat(paste("-pdf saved at", output, "\n"))
}


main <- function(args){
  # load data
  df_cln <- load_table(args$cln_table)
  df_fus <- load_table(args$fus_table)
  df_sam <- load_table(args$sam_table)

  # select samples
  df_cln_rna <- df_cln %>% filter(!is.na(Sample_Id_RNA_T))
  df_fus <- df_fus %>% filter(Sample_Id %in% df_cln_rna$Sample_Id_RNA_T)
  df_fus <- df_fus %>% filter(Sample_Id %in% df_sam$Sample_Id)

  # draw without breakpoints
  df_fus_wo_breakpoint <- df_fus %>% select(Sample_Id, Fusion_Id, FILTER) %>% distinct()
  draw_upset_plot_fusion(df_fus=df_fus_wo_breakpoint, output=args$output_wo_breakpoints, field="FILTER")
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregated fusion calls for one cohort across different callers.')
  parser$add_argument("--cln_table", type="character",
                      help="Path to table of aggregated and annotated prism fusion calls.",
                      default="../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
  parser$add_argument("--fus_table", type="character",
                      help="Path to table of aggregated and annotated prism fusion calls.",
                      default="../../../data/prism/rna/fusions/prism_filters.tsv.gz")
  parser$add_argument("--sam_table", type="character",
                      help="Path to table of samples.",
                      default="../../../results/fusions_analysis/selection/selection_samples_prism.tsv")
  parser$add_argument("--output_wo_breakpoints", type="character",
                      help="Upset plot of filtered fusions without breakpoint.",
                      default="../../../results/fusions_analysis/other_plots/upset_callers_wt_breakpoints_prism.pdf")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
