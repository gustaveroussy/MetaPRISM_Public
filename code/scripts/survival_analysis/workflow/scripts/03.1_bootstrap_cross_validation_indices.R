# @created: 13 Apr 22
# @modified: 05 May 22
# @authors: Yoann Pradat
#
#    CentraleSupelec
#    MICS laboratory
#    9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
#
#    Institut Gustave Roussy
#    Prism Center
#    114 rue Edouard Vaillant, Villejuif, 94800 France
#
# Generate a table with lists of indices from 0 to the maximum of observations - 1 for performing bootstrap and nested
# cross-validations in the following format.
# 
# Run       |  Index    |  Split
# 
# Boot_1    |    0      |  Train
# Boot_1    |    8      |  Test
# Boot_1    |    10     |  Train
#   ...
# Boot_100  |   164     |  Test
# 
# Generate another table with list of indices from 0 to the maximum of observations - 1 for performing repeated
# cross-validations and nested cross-validations.
# 
# 
# Run       |  Index    |  Split
# 
# CV_Train  |    0      |  1
# CV_1      |    8      |  2
# CV_1      |    10     |  2
#   ...
# CV_200    |   164     |  5
# 
# 
# Bootstrap samples and cross-validation are selected so that the proportion of each strata of the outcome is stable
# across samples/splits.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

# functions ============================================================================================================

generate_boot_indices <- function(indices_all, n_boot){
  df_indices <- tibble()

  for (b in 1:n_boot){
    train_indices <- sample(indices_all, size=length(indices_all), replace=T)
    test_indices <- setdiff(indices_all, train_indices)
    df_indices_train <- tibble(Run=paste0("Boot_", b), Index=train_indices, Split="Train")
    df_indices_test <- tibble(Run=paste0("Boot_", b), Index=test_indices, Split="Test")
    df_indices_b <- bind_rows(df_indices_train, df_indices_test)
    df_indices <- bind_rows(df_indices, df_indices_b)
  }

  df_indices
}


generate_xval_indices <- function(indices_all, k_cv, run){
  df_indices <- tibble()
  indices_left <- indices_all
  test_size <- round(length(indices_all)/k_cv)

  for (k in 1:k_cv){
    if (k==k_cv){
      indices_k <- indices_left
    } else {
      indices_k <- sample(indices_left, size=test_size, replace=F)
    }

    indices_left <- setdiff(indices_left, indices_k)
    df_indices <- bind_rows(df_indices, tibble(Run=run, Index=indices_k, Split=k))
  }

  df_indices
}


save_table <- function(table, filepath){
  # save
  if (grepl(".gz$", filepath)){
    write.table(table, gsub(".gz$", "", filepath), sep="\t", row.names=F, quote=F)
    system(paste("gzip", gsub(".gz$", "", filepath)))
  } else {
    write.table(table, filepath, sep="\t", row.names=F, quote=F)
  }
  cat("-table saved at", filepath, "\n")
}


main <- function(args){
  # identify path to data table
  if ("data.imputed_1.tsv.gz" %in% list.files(args$dat_folder)){
    dat_table <- file.path(args$dat_folder, "data.imputed_1.tsv.gz")
  } else if ("data.complete_cases.tsv.gz" %in% list.files(args$dat_folder)){
    dat_table <- file.path(args$dat_folder, "data.complete_cases.tsv.gz")
  } else {
    stop(paste("-ERROR: the folder", args$dat_folder, "doest not contain any of data.imputed_1.tsv.gz or",
               "data.complete_cases.tsv.gz"))
  }

  # read  data
  df_dat <- read_tsv(dat_table, show_col_types=F, progress=F)
  n_row <- nrow(df_dat)

  # for reproducibility
  set.seed(args$seed)

  # generate bootstrap indices =========================================================================================
  df_indices_b <- generate_boot_indices(indices_all=1:n_row, n_boot=args$n_boot)
  df_indices_b_full <- generate_xval_indices(indices_all=1:n_row, k_cv=args$k_cv, run="Full")
  df_indices_b_full$Split <- as.character(df_indices_b_full$Split)
  df_indices_b <- bind_rows(df_indices_b_full, df_indices_b)

  # save
  save_table(df_indices_b, args$output_boot)

  # generate cross-validation indices ==================================================================================
  df_indices_x <- tibble()

  # repeated cross-validation
  for (x in 0:args$n_xval){
    if (x==0){
      run <- "Full"
    } else {
      run <- paste0("CV_", x)
    }

    df_indices <- generate_xval_indices(indices_all=1:n_row, k_cv=args$k_cv, run=run)
    
    # nested cross-validations
    for (split in unique(df_indices$Split)){
      df_indices_split <- df_indices %>% filter(Split!=split) %>% select(Run, Index, Split)
      indices_split <- df_indices_split %>% pull(var="Index")
      df_indices_split_nested <- generate_xval_indices(indices_all=indices_split, k_cv=args$k_cv, run=run)
      df_indices_split_nested <- df_indices_split_nested %>% rename(!!paste0("Split_Nested", "_", split):=Split)
      df_indices_split <- left_join(df_indices_split, df_indices_split_nested, by=c("Run", "Index"))
      df_indices <- left_join(df_indices, df_indices_split, by=c("Run", "Index", "Split"))
    }

    df_indices_x <- bind_rows(df_indices_x, df_indices)
  }

  # save
  save_table(df_indices_x, args$output_xval)
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder <- "../../../results/survival_analysis/data/sub_features/BRCA_dna_and_rna/cln_biol_0/processed"

  parser <- ArgumentParser(description='Prepare tables with indices for bootstrapping or cross-validation.')
  parser$add_argument('--dat_folder', type="character", default=default_folder,
                      help='Folder in which the data table will be searched for$')
  parser$add_argument('--n_boot', type="integer", default=100, help='Number of bootstrap samples.')
  parser$add_argument('--n_xval', type="integer", default=200, help='Number of cross-validation repeats.')
  parser$add_argument('--k_cv', type="integer", default=5, help='Number of folds for nested cross-validation.')
  parser$add_argument('--seed', type="integer", default=1995, help='Seed for the random number generator.')
  parser$add_argument('--output_boot', type="character", help='Path to where table of bootstrap indices will be saved.')
  parser$add_argument('--output_xval', type="character",
                      help='Path to where table of cross-validation indices will be saved.')
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
