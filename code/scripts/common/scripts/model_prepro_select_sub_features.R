# @created: 05 Nov 21
# @modified: 05 Aug 22
# @authors: Yoann Pradat
#
# Considering the large table produced in the script model_prepro_gather_all_features.R, prepare sub tables with 
# selected columns and selected rows.
#
#  - Columns selection is driven by the args$features parameter
#  - Rows selection is driven by the args$features and args$samples parameters. Features intervene here as in case 
#  DNA-seq and/or RNA-seq features are used, only patients for which DNA and/or RNA was sequenced may be used.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(yaml))

source("../common/functions/model_utils.R")

Sys.setenv("VROOM_CONNECTION_SIZE"=524288)

# functions ============================================================================================================

get_covs_list <- function(covs, features, features_all, df_cov){
  for (cov in features){
    if (cov %in% names(features_all)){
      covs <- get_covs_list(covs, features_all[[cov]], features_all, df_cov)
    } else {
      covs_more <- df_cov %>% filter(grepl(cov, df_cov$Class)) %>% pull(Covariate)
      if (length(covs_more) == 0){
        covs_more <- df_cov %>% filter(df_cov$Covariate==cov) %>% pull(Covariate)
      }
      covs <- c(covs, covs_more)
    }
  }

  covs
}


add_features_indicators <- function(df_cov, features_all){
  for (sub_name in names(features_all)){
    covs <- get_covs_list(covs=c(), features=sub_name, features_all, df_cov)
    df_cov <- df_cov %>% mutate(!!sub_name:=ifelse(Covariate %in% covs, 1, 0))
  }
  
  df_cov
}


select_rows_samples <- function(df_dat, df_cnt, df_sam, samples){
  df_dat_sub <- NULL

  # select tumor types to keep from value of "samples"
  tt_keep <- NULL

  for (pattern in c("all_cohort", "dna_cohort", "rna_cohort")){
    if (grepl(paste0("^", pattern), samples)){
      col_usage <- paste0("Use", "_", pattern)
      tt_keep <- df_cnt %>% filter(.data[[col_usage]]==1) %>% pull(Tumor_Type)
    }
  }

  if (is.null(tt_keep)){
    tt_all <- df_dat %>% pull(Project_TCGA_More)
    for (tt in tt_all){
      if (grepl(paste0("^", tt), samples)){
        tt_keep <- tt
      }
    }
  }

  if (is.null(tt_keep)){
    stop(paste("-error! could not determine which tumor types should be kept from the value", samples))
  }

  df_dat_sub <- df_dat %>% filter(Project_TCGA_More %in% tt_keep)
  df_sam_sub <- df_sam %>% filter(Tumor_Type %in% tt_keep)

  # select possibly only samples dna and/or rna 
  regexes <- list(dna="DNA_T", rna="RNA_T", dna_and_rna="DNA_T\\|RNA_T", dna_or_rna="DNA_T|RNA_T")

  for (pattern in c("dna_and_rna", "dna_or_rna", "dna", "rna")){
    regex <- regexes[[pattern]]
    if (grepl(paste0(pattern, "$"), samples)){
      col_sam <- paste("Use", pattern, "cohort", sep="_")
      df_sam_sub <- df_sam_sub %>% filter(.data[[col_sam]]==1)
      df_dat_sub <- df_dat_sub %>% filter(Subject_Id %in% df_sam_sub$Subject_Id)
      # df_dat_sub <- df_dat_sub %>% filter(grepl(regex, Sample_Type, perl=T))
      break
    }
  }

  df_dat_sub
}


select_rows_features <- function(df_dat, df_cov_sub){
  df_dat_sub <- df_dat

  covs_classes <- df_cov_sub %>% pull(Class) %>% unique()
  dna_used <- any(grepl("^DNA", covs_classes))
  rna_used <- any(grepl("^RNA", covs_classes))

  if (dna_used){
    nrow_b <- nrow(df_dat_sub)
    df_dat_sub <- df_dat_sub %>% filter(grepl("DNA_N\\|DNA_T", Sample_Type))
    cat(paste("-selected", paste0(nrow(df_dat_sub), "/", nrow_b), "samples with DNA data\n"))
  }
  if (rna_used){
    nrow_b <- nrow(df_dat_sub)
    df_dat_sub <- df_dat_sub %>% filter(grepl("RNA_T", Sample_Type))
    cat(paste("-selected", paste0(nrow(df_dat_sub), "/", nrow_b), "samples with RNA data\n"))
  }

  # special treatment for mutational signatures as samples with too low mutational burden were not analyzed
  mut_sig_used <- any(grepl("Mutational_Signatures", covs_classes))
  if (mut_sig_used){
    nrow_b <- nrow(df_dat_sub)
    covs_sig <- df_cov_sub %>% filter(grepl("Mutational_Signatures", Class)) %>% pull(Covariate)
    df_dat_sub <- df_dat_sub %>% filter(!is.na(.data[[covs_sig[1]]]))
    cat(paste("-selected", paste0(nrow(df_dat_sub), "/", nrow_b), "samples with DNA mutational signatures data\n"))
  }


  # special treatment for Bagaev tme as some tumor types were not analyzed by Bagaev.
  tme_bagaev_used <- any(grepl("TME_Bagaev", covs_classes))
  if (tme_bagaev_used){
    nrow_b <- nrow(df_dat_sub)
    covs_tme <- df_cov_sub %>% filter(grepl("TME_Bagaev", Class)) %>% pull(Covariate)
    df_dat_sub <- df_dat_sub %>% filter(!is.na(.data[[covs_tme[1]]]))
    cat(paste("-selected", paste0(nrow(df_dat_sub), "/", nrow_b), "samples with RNA TME Bagaev data\n"))
  }

  # special treatment for fusions as some fusion calling algorithms failed on some tumor samples
  fus_used <- any(grepl("RNA-Fusion", covs_classes))
  if (fus_used){
    nrow_b <- nrow(df_dat_sub)
    covs_fus <- df_cov_sub %>% filter(grepl("RNA-Fusion", Class)) %>% pull(Covariate)
    df_dat_sub <- df_dat_sub %>% filter(!is.na(.data[[covs_fus[1]]]))
    cat(paste("-selected", paste0(nrow(df_dat_sub), "/", nrow_b), "samples with RNA fusion data\n"))
  }

  df_dat_sub
}


main <- function(args){
  # load data and covariates tables
  df_dat <- load_table(args$input_dat)
  df_cov <- load_table(args$input_cov)
  df_cnt <- load_table(args$input_cnt)
  df_sam <- load_table(args$input_sam)
  names_covs <- names(df_cov)
  
  cat(paste("-this script will perform columns and rows selection for", args$features, "and", args$samples, "\n"))

  # read features from yaml
  features_all <- read_yaml(args$config_yaml)[[args$config_section]]$features

  # build class and add sub features indicators
  df_cov <- unite_class_levels(df_cov)
  df_cov <- add_features_indicators(df_cov, features_all)

  # outcome covs
  covs_outcome <- df_cov %>% filter(grepl("Outcome", Class) & !grepl("Other_Outcome", Class)) %>% pull(Covariate)

  # select covs
  df_cov_sub <- df_cov %>% filter(.data[[args$features]]==1)
  covs <- df_cov_sub %>% pull(Covariate)
  covs <- c("Subject_Id", covs_outcome, covs)
  cat(paste("-selected", length(covs), "covariates\n"))

  # select rows
  df_dat_sub <- select_rows_samples(df_dat, df_cnt, df_sam, args$samples)
  df_dat_sub <- select_rows_features(df_dat_sub, df_cov_sub)

  # add cov outcomes to cov_sub
  df_cov_out <- df_cov %>% filter(Covariate %in% covs_outcome)
  df_cov_sub <- bind_rows(df_cov_out, df_cov_sub)

  # keep only names originally present in covariates table 
  df_dat_sub <- df_dat_sub[,covs]
  df_cov_sub <- df_cov_sub %>% select(any_of(names_covs))

  # save
  dir.create(gsub(basename(args$output_dat), "", args$output_dat), showWarnings=F, recursive=T)
  save_table(df_dat_sub, args$output_dat, sep="\t")
  save_table(df_cov_sub, args$output_cov, sep="\t")
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Prepare tables with sub features selected.')
  parser$add_argument("--input_dat", type="character", help="Path to input data table.",
                      default="../../../results/survival_analysis/data_tcga/all_features/data.tsv.gz")
  parser$add_argument("--input_cov", type="character", help="Path to input covariates table.",
                      default="../../../results/survival_analysis/data_tcga/all_features/covs.tsv.gz")
  parser$add_argument("--input_sam", type="character", help="Path to table of samples selection.",
                      default="../../../results/survival_analysis/selection/selection_samples_tcga.tsv")
  parser$add_argument("--input_cnt", type="character", help="Path to input table of tumor type counts.",
                      default="../../../results/survival_analysis/selection/selection_tumor_types.tsv")
  parser$add_argument("--config_yaml", type="character", default="config/config.yaml",
                      help="Path to the config.yaml file where subselections of features are defined.")
  parser$add_argument("--config_section", type="character", default="models",
                      help="Name of the section containing parameters for the models.")
  parser$add_argument("--features", type="character", help="Name of the selection of features.",
                      default="cln_dna_1_rna_1")
  parser$add_argument("--samples", type="character", help="Name of the selection of samples",
                      default="BRCA_dna_and_rna")
  parser$add_argument("--output_dat", type="character", help="Path to output data table")
  parser$add_argument("--output_cov", type="character", help="Path to output covariates table")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
