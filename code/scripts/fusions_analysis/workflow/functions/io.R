suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

load_samples_for_analysis <- function(cohorts, samples, analysis=NULL){
  dfs_sam <- setNames(lapply(samples, read_tsv, col_types=cols()), cohorts)

  for (cohort in cohorts){
    df_sam <- dfs_sam[[cohort]]
    if (!is.null(analysis)){
      df_sam <- df_sam %>% filter(.data[[paste0("Use_For_", analysis)]] == 1)
    }
    dfs_sam[[cohort]] <- df_sam
  }

  dfs_sam
}


load_fusions_for_analysis <- function(cohorts, fusions, samples=NULL, analysis=NULL){
  dfs_fus <- setNames(lapply(fusions, read_tsv, col_types=cols()), cohorts)

  if (!is.null(samples)){
    if (is.character(samples[[1]])){
      dfs_sam <- load_samples_for_analysis(cohorts, samples, analysis)
    } else {
      dfs_sam <- samples
    }
  }

  for (cohort in cohorts){
    df_fus <- dfs_fus[[cohort]]

    if (!is.null(analysis)){
      df_fus <- df_fus %>% filter(.data[[paste0("Use_For_", analysis)]] == 1)
    }
    
    if (!is.null(samples)){
      df_sam <- dfs_sam[[cohort]]
      df_fus <- df_fus %>% filter(Sample_Id %in% df_sam$Sample_Id)
    }

    dfs_fus[[cohort]] <- df_fus
  }

  dfs_fus
}
