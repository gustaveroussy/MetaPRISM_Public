# @created: 11 Feb 22
# @modified: 11 Feb 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

load_table_samples <- function(algo, cohort){
  filepath <- paste0("../../../data/", cohort, "/rna/", algo, "/sample_list.tsv")
  df_sam <- read_tsv(filepath, progress=F, show_col_types=F)
  df_sam["Algo"] <- algo
  df_sam[,c("Sample_Id", "Algo")]
}


aggregate_callers_per_id <- function(df_fus, col_id){
  cols_call <- colnames(df_fus)[grepl("_Call$", colnames(df_fus))]
  dfs_fus_call <- list()
  for (col_call in cols_call){
    dfs_fus_call[[col_call]] <- df_fus %>% filter(.data[[col_call]]==1) %>% select(all_of(c("Sample_Id", col_id))) %>%
      distinct() %>% mutate(Caller=col_call)
  }
  df_fus_call <- bind_rows(dfs_fus_call)

  df_fus_call %>% group_by_at(all_of(c("Sample_Id", col_id))) %>%
    summarize(Caller=paste0(Caller, collapse="|"), .groups="drop")
}


extract_commonly_analyzed_samples <- function(df_sam_tcga, df_sam_tcga_validation){
  algos_tcga <- sort(unique(df_sam_tcga$Algo))
  algos_tcga_validation <- sort(unique(df_sam_tcga_validation$Algo))

  df_sam_tcga_agg <- df_sam_tcga %>% arrange(Algo) %>% group_by(Sample_Id) %>%
    summarize(Algo=paste0(Algo, collapse="|"), .groups="drop")
  df_sam_tcga_validation_agg <- df_sam_tcga_validation %>% arrange(Algo) %>% group_by(Sample_Id) %>%
    summarize(Algo=paste0(Algo, collapse="|"), .groups="drop")

  sam_tcga_all <- df_sam_tcga_agg %>% filter(Algo==paste0(algos_tcga, collapse="|")) %>% pull(Sample_Id)
  sam_tcga_validation <- df_sam_tcga_validation_agg %>% filter(Algo==paste0(algos_tcga_validation, collapse="|")) %>%
    pull(Sample_Id)

  intersect(sam_tcga_all, sam_tcga_validation)
}


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


draw_upset_plot_sample <- function(df_sam, output){
  df_sam_agg <- df_sam %>% group_by(Sample_Id) %>%
    summarize(Algo=paste0(Algo, collapse="|"), .groups="drop")

  df_sam_indices <- build_indices(df=df_sam_agg, field="Algo", unique_values=unique(df_sam$Algo))
  m_upset <- make_upset_m(df_sam_indices, field_set="Set", field_identifier="Index")
  pdf(file=output, width=8, height=4)
  draw_upset_plot(m_upset, width_set_size=2, height_top_annot=2)
  dev.off()
  cat(paste("-pdf saved at", output, "\n"))
}



draw_upset_plot_fusion <- function(df_fus, col_id,  output){
  cols_call <- colnames(df_fus)[grepl("_Call$", colnames(df_fus))]
  df_fus_agg <- aggregate_callers_per_id(df_fus, col_id=col_id)
  df_fus_indices <- build_indices(df=df_fus_agg, field="Caller", unique_values=cols_call)
  m_upset <- make_upset_m(df_fus_indices, field_set="Set", field_identifier="Index")
  pdf(file=output, width=12, height=6)
  draw_upset_plot(m_upset, width_set_size=3, height_top_annot=3)
  dev.off()
  cat(paste("-pdf saved at", output, "\n"))
}


main <- function(args){
  # load samples
  df_sam_tcga <- bind_rows(lapply(args$algos_tcga, load_table_samples, cohort="tcga"))
  df_sam_tcga_validation <- bind_rows(lapply(args$algos_tcga_validation, load_table_samples, cohort="tcga_validation"))

  # upset plot samples
  if (args$cohort=="tcga"){
    draw_upset_plot_sample(df_sam=df_sam_tcga, output=args$output_samples)
  } else if (args$cohort=="tcga_validation"){
    draw_upset_plot_sample(df_sam=df_sam_tcga_validation, output=args$output_samples)
  }

  # extract commonly analyzed samples
  sam_com <- extract_commonly_analyzed_samples(df_sam_tcga=df_sam_tcga,
                                               df_sam_tcga_validation=df_sam_tcga_validation)

  # load table
  df_fus <- read_tsv(args$fusions, show_col_types=F, progress=F)

  # reset FFPM >= 0.1 for starfusion
  if ("SF_Call" %in% colnames(df_fus)){
    df_fus <- df_fus %>% mutate(SF_Call=ifelse(SF_Call==0,0,ifelse(SF_FFPM>=0.1, 1, 0))) 
  }

  # add ids with breakpoint
  df_fus <- df_fus %>% 
    mutate(Fusion_Id_With_Breakpoint=paste(Fusion_Id, "_", paste0(Breakpoint_1, "--", Breakpoint_2)))

  # split between all samples and commonly analyzed samples
  df_fus_all <- df_fus
  df_fus_com <- df_fus %>% filter(Sample_Id %in% sam_com)

  # upset plot per Fusion_Id and Fusion_Id_With_Breakpoint
  draw_upset_plot_fusion(df_fus=df_fus_all, col_id="Fusion_Id",
                         output=args$output_wo_breakpoints_all)
  draw_upset_plot_fusion(df_fus=df_fus_com, col_id="Fusion_Id",
                         output=args$output_wo_breakpoints_com)
  draw_upset_plot_fusion(df_fus=df_fus_all, col_id="Fusion_Id_With_Breakpoint",
                         output=args$output_w_breakpoints_all)
  draw_upset_plot_fusion(df_fus=df_fus_com, col_id="Fusion_Id_With_Breakpoint",
                         output=args$output_w_breakpoints_com)
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregated fusion calls for one cohort across different callers.')
  parser$add_argument("--fusions", type="character",
                      help="Path to table of aggregated and annotated tcga fusion calls.",
                      default="../../../data/tcga/rna/fusions/tcga_annotated.tsv.gz")
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="tcga")
  parser$add_argument("--algos_tcga", type="character", nargs="+", help="Algos analyzed for tcga data",
                      default=c("deepest_pnas_2019", "prada_nar_2018", "starfusion_cell_2018"))
  parser$add_argument("--algos_tcga_validation", type="character", nargs="+",
                      help="Algos analyzed for tcga validation data",
                      default=c("arriba", "ericscript", "fusioncatcher", "pizzly", "squid", "starfusion"))
  parser$add_argument("--output_samples", type="character", help="Upset plot samples tcga",
          default="../../../results/fusions_analysis/tcga_validation/upset_samples_tcga.pdf")
  parser$add_argument("--output_wo_breakpoints_all", type="character", help="Upset plot all fusions.",
          default="../../../results/fusions_analysis/tcga_validation/upset_callers_wo_breakpoints_all_tcga.pdf")
  parser$add_argument("--output_wo_breakpoints_com", type="character",
                      help="Upset plot fusions of commonly analyzed samples.",
          default="../../../results/fusions_analysis/tcga_validation/upset_callers_wo_breakpoints_com_tcga.pdf")
  parser$add_argument("--output_w_breakpoints_all", type="character", help="Upset plot all fusions.",
          default="../../../results/fusions_analysis/tcga_validation/upset_callers_w_breakpoints_all_tcga.pdf")
  parser$add_argument("--output_w_breakpoints_com", type="character",
                      help="Upset plot fusions of commonly analyzed samples.",
          default="../../../results/fusions_analysis/tcga_validation/upset_callers_w_breakpoints_com_tcga.pdf")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
