# @created: 20 Oct 21
# @modified: 06 Dec 22
# @authors: Yoann Pradat
#
# Describe the distribution of drugs before biopsy according to tumor types.

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

add_drug_class <- function(df_drug_table, max_level=NULL){
  levels = sort(colnames(df_drug_table)[grepl("Class_Lvl", colnames(df_drug_table))])
  if (!is.null(max_level)){
    levels <- levels[1:max_level]
  } 
  df_drug_table <- df_drug_table %>% unite("Class", all_of(levels), sep=" - ", remove=F) %>%
    mutate(Class=gsub(" - NA", "", Class))
  df_drug_table
}


main <- function(args){
  # load clinical data
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  col_row <- args$col_row

  # load drug table
  df_drug_table <- load_table(args$drug_table)
  df_drug_table <- add_drug_class(df_drug_table)

  # group rare tumor types
  rare_type <- "Rare subtypes"
  unknown_type <- "Unknown primary"
  df_cln <- df_cln %>% mutate(Tumor_Type=ifelse(grepl("Not_TCGA", Tumor_Type), rare_type, Tumor_Type))
  df_cln <- df_cln %>% mutate(Tumor_Type=ifelse(grepl("Unknown", Tumor_Type), unknown_type, Tumor_Type))

  # list tumor types represented by at least 10 tumors
  tumor_types <- df_cln %>% group_by(Tumor_Type) %>% summarize(n=n()) %>% filter(n>=args$min_count_cols) %>% 
    arrange(desc(n)) %>% pull(Tumor_Type)
  tumor_types <- c(setdiff(tumor_types, c(rare_type, unknown_type)), c(rare_type, unknown_type))
  
  # filter out patients with no treatments data
  df_cln <- df_cln %>% filter(!is.na(.data[[col_row]]))
  # mask_na <- is.na(df_cln[["Drugs_Before_Biopsy"]])
  # mask_com <- df_cln[["Comment_Drugs_Before_Biopsy"]] %in% c("Pas de dossier", "Pas de traitement prébiopsie")
  # mask_rm <- mask_na & mask_com
  # df_cln <- df_cln %>% filter(!mask_rm)
  # df_cln <- df_cln %>% filter(.data[[col_row]]!="None")

  # plot all tumor types
  df_count <- df_cln %>% separate_rows(.data[[col_row]], sep="\\|") %>% 
    select(Subject_Id, Tumor_Type, all_of(col_row)) %>% filter(!is.na(.data[[col_row]]))
  df_count <- df_count %>% group_by_at(all_of(col_row)) %>% filter(n()>=args$min_count_rows)
  df_count <- dcast(df_count, formula = Subject_Id + Tumor_Type ~ df_count[[col_row]], fun.aggregate = length)

  plot_data <- get_plot_data(df_count, dcolor_quantile=0.5)
  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[[1]],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors that received the drug",
                   dcolor_title_legend="Median number of drugs received")

  # plot selected tumor types
  df_count <- df_cln %>% filter(Tumor_Type %in% tumor_types) %>%
    separate_rows(.data[[col_row]], sep="\\|") %>%
    select(Subject_Id, Tumor_Type, all_of(col_row)) %>% filter(!is.na(.data[[col_row]]))
  
  rows_min_count <- df_count %>% group_by_at(all_of(col_row)) %>% filter(n()>=args$min_count_rows) %>% pull(col_row)
  rows_min_count <- unique(rows_min_count)
  rows_min_frac <- df_count %>% group_by_at(all_of(c(col_row, "Tumor_Type"))) %>%
    mutate(N_Drug=n()) %>% group_by(Tumor_Type) %>% mutate(N_Tumor=length(unique(Subject_Id))) %>%
    mutate(Frac=N_Drug/N_Tumor) %>% filter(Frac>=args$min_frac_rows) %>% pull(col_row)
  rows_min_frac <- unique(rows_min_frac)
  rows_keep <- union(rows_min_count, rows_min_frac)
  rows_keep <- setdiff(rows_keep, "None")

  # order rows by class
  if (args$col_row!="Classes_Before_Biopsy"){
    classes_keep <- df_drug_table$Class[match(rows_keep, df_drug_table$DCI)]
    df_rows_keep <- tibble(Drug=rows_keep, Class=classes_keep)

    classes_simplify <- list("Antiangiogenic - VEGFA/VEGFB/VEGFC/VEGFD"="Antiangio.",
                             "Chemo_Alkylating"="Chemo. Alk.",
                             "Chemo_Non_Alkylating"="Chemo. Non-Alk.",
                             "Hormonotherapy - Antiandrogen"="Hormono. AR",
                             "Hormonotherapy - Antioestrogen - Aromatase inhibitor - Non-steroidal"="Hormono. Oes",
                             "Hormonotherapy - Antioestrogen - Aromatase inhibitor - Steroidal"="Hormono. Oes",
                             "Hormonotherapy - Antioestrogen - Selective estrogen receptor degrader - Steroidal"="Hormono. Oes",
                             "Hormonotherapy - Antioestrogen - Selective estrogen receptor modulator"="Hormono. Oes",
                             "Immunotherapy - Checkpoint_inhibitor - PD-L1"="Immuno.",
                             "Immunotherapy - Checkpoint_inhibitor - PD1"="Immuno.",
                             "mTORi"="i. mTOR",
                             "Platinum_salts"="A. Platinum Salts",
                             "Targeted_Therapy - ALK/ROS1"="i. ALK",
                             "Targeted_Therapy - EGFR"="i. EGFR",
                             "Targeted_Therapy - FLT1/KDR/PDGFRA/PDGFRB/KIT/FLT3/CSF1R/RET"="i. RTKs",
                             "Targeted_Therapy - HER2"="i. HER2",
                             "Targeted_Therapy - KDR"="i. KDR",
                             "Targeted_Therapy - RAF"="i. RAF")
    df_rows_keep$Class_Simple <- sapply(match(df_rows_keep$Class, names(classes_simplify)),
                                        function(i) classes_simplify[[i]])
    df_rows_keep <- df_rows_keep %>% arrange(Class_Simple, Drug)
    rows_keep <- df_rows_keep$Drug
  } else {
    rows_keep <- sort(rows_keep)
  }

  # shorten some drug names
  drugs_shorten <- list("DOXORUBICINE LIPOSOMALE"="DOXORUBICINE LIPO.")
  for (i in 1:length(drugs_shorten)){
    drug_cur <- names(drugs_shorten)[[i]]
    drug_new <-  drugs_shorten[[i]]
    df_count[[col_row]] <- gsub(drug_cur, drug_new, df_count[[col_row]])
    rows_keep <- gsub(drug_cur, drug_new, rows_keep)
  }

  df_count <- df_count %>% filter(.data[[col_row]] %in% rows_keep)
  df_count <- dcast(df_count, formula = Subject_Id + Tumor_Type ~ df_count[[col_row]], fun.aggregate = length)
  df_tt_count <- df_cln %>% group_by(Tumor_Type) %>% summarize(n=n())
  df_count <- left_join(df_count, df_tt_count, by="Tumor_Type")
  df_count$Tumor_Type <- factor(df_count$Tumor_Type, levels=tumor_types)
  df_count <- df_count[,c("Subject_Id", "Tumor_Type", rows_keep, "n")]

  plot_data <- get_plot_data(df_count, dcolor_quantile=0.5)
  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[[2]],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors that received the drug",
                   dcolor_title_legend="Median number of drugs received")

  if (!is.null(args$output_paper)){
    draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$output_paper,
                     dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                     rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                     dscale_title_legend="Proportion of tumors that received the drug",
                     dcolor_title_legend="Median number of drugs received")
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe distribution of drugs before biopsy per tumor type.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument('--drug_table', type="character", help='Path to the drug table.',
                      default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
  parser$add_argument('--min_count_cols', type="integer", default=10,
                      help='Each column (i.e tumor type) with at least --min_count_cols tumors will be in the plot.')
  parser$add_argument('--min_count_rows', type="integer", default=20,
                      help='Each drug received by at least --min_count_rows patients overall will be in the plot.')
  parser$add_argument('--min_frac_rows', type="double", default=0.5,
                      help='Each drug received by at least --min_frac_rows of any tumor type will be in the plot.')
  parser$add_argument('--col_row', type="character", default="Drugs_Before_Biopsy",
                      help='Choose between Drugs_Before_Biopsy and Classes_Before_Biopsy.')
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plots.",
    default=c("../../../results/data_overview/treatments/heatmaps/table_extra_drugs_Drugs_Before_Biopsy_all_tt.pdf",
              "../../../results/data_overview/treatments/heatmaps/table_extra_drugs_Drugs_Before_Biopsy_sel_tt.pdf"))
  parser$add_argument("--output_paper", type="character", help="Path to output plots.",
                      default=NULL)
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
