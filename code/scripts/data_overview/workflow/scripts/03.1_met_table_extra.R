# @created: 16 Oct 21
# @modified: 29 Oct 21
# @authors: Yoann Pradat
#
# Describe the distribution of metastatic sites according to tumor types.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tableExtra))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

source("workflow/functions/table_extra.R")

# functions ============================================================================================================

main <- function(args){
  # load clinical data
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)

  # plot all tumor types
  df_count <- df_cln %>% separate_rows(Metastatic_Sites, Metastatic_Subsites, sep="\\|") %>%
    select(Subject_Id, Tumor_Type, Metastatic_Sites) %>% filter(!is.na(Metastatic_Sites))
  df_count <- dcast(df_count, formula = Subject_Id + Tumor_Type ~ Metastatic_Sites, fun.aggregate = length)

  plot_data <- get_plot_data(df_count, dcolor_quantile=0.5)
  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[[1]],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more,
                   dscale_title_legend="Proportion of tumors metastatic at the site",
                   dcolor_title_legend="Median number of metastatic sites")

  # plot selected tumor types
  tumor_types <- df_cln %>% group_by(Tumor_Type) %>% summarize(n=n()) %>% filter(n>=args$min_count) %>%
    filter(!Tumor_Type %in% c("MISC - Not_TCGA")) %>% pull(Tumor_Type)
  df_count <- df_cln %>% filter(Tumor_Type %in% tumor_types) %>%
    separate_rows(Metastatic_Sites, Metastatic_Subsites, sep="\\|") %>%
    select(Subject_Id, Tumor_Type, Metastatic_Sites) %>% filter(!is.na(Metastatic_Sites))
  df_count <- dcast(df_count, formula = Subject_Id + Tumor_Type ~ Metastatic_Sites, fun.aggregate = length)

  plot_data <- get_plot_data(df_count, dcolor_quantile=0.5)
  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[[2]],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more,
                   dscale_title_legend="Proportion of tumors metastatic at the site",
                   dcolor_title_legend="Median number of metastatic sites")

  # plot selected tumor types + group into Rare subtypes
  # group rare tumor types
  rare_type <- "Rare subtypes"
  unknown_type <- "Unknown primary"
  df_cln <- df_cln %>% mutate(Tumor_Type=ifelse(grepl("Not_TCGA", Tumor_Type), rare_type, Tumor_Type))
  df_cln <- df_cln %>% mutate(Tumor_Type=ifelse(grepl("Unknown", Tumor_Type), unknown_type, Tumor_Type))
  tumor_types <- df_cln %>% group_by(Tumor_Type) %>% summarize(n=n()) %>% filter(n>=args$min_count) %>%
    arrange(desc(n)) %>% pull(Tumor_Type)
  df_count <- df_cln %>% filter(Tumor_Type %in% tumor_types) %>%
    separate_rows(Metastatic_Sites, Metastatic_Subsites, sep="\\|") %>%
    select(Subject_Id, Tumor_Type, Metastatic_Sites) %>% filter(!is.na(Metastatic_Sites))
  df_count <- dcast(df_count, formula = Subject_Id + Tumor_Type ~ Metastatic_Sites, fun.aggregate = length)

  # order tumor types
  tumor_types <- c(setdiff(tumor_types, c(rare_type, unknown_type)), c(rare_type, unknown_type))
  df_count$Tumor_Type <- factor(df_count$Tumor_Type, levels=tumor_types)

  plot_data <- get_plot_data(df_count, dcolor_quantile=0.5)
  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$output_paper,
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more,
                   dscale_title_legend="Proportion of tumors metastatic at the site",
                   dcolor_title_legend="Median number of metastatic sites")

}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe distribution of metastatic sites per tumor type.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument('--min_count', type="integer", default=10,
                      help="Minimum number of tumors for a tumor type to be in restricted plots.")
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plots.",
    default=c("../../../results/data_overview/metastatic_sites/table_extra_met_sites_all_tt.pdf",
              "../../../results/data_overview/metastatic_sites/table_extra_met_sites_sel_tt.pdf"))
  parser$add_argument("--output_paper", type="character", help="Path to output paper plots.",
                      default="../../../results/figures_paper/FS2.pdf")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
