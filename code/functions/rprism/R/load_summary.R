# @created: 30 Aug 21
# @modified: 30 Aug 21
# @author: Yoann Pradat
# 
# Exported functions
#   load_summary_rna_fus
#   load_summary_rna_gex
#   load_summary_wes_mut

#' Load summary of rna-seq gene fusion tables
#' 
#' The summary table links Sample_Id to Subject_Id. This allows easier connection
#' to other data tables using either of Sample_Id or Subject_Id.
#'
#' @inheritParams load_rna_fus
#' @return a data.frame with summary
#' @export
load_summary_rna_fus <- function(study, mode=NULL){
  load_from_data(get_filepath_summary_rna_fus(study, mode), progress=F, show_col_types=F)
}


#' Load summary of RNA tables
#' 
#' The summary table links original names in the RNA table to Sample_Id and Subject_Id. This allows easier connection
#' to other data tables using either of Sample_Id or Subject_Id.
#'
#' @inheritParams load_rna_gex
#' @return a data.frame with summary
#' @export
load_summary_rna_gex <- function(study, level, metric, test_mode=F){
  load_from_data(get_filepath_summary_rna_gex(study, level, metric, test_mode), progress=F, show_col_types=F)
}

#' Load summary of MAF tables
#' 
#' The summary table links original names in the WES table to Sample_Id and Subject_Id. This allows easier connection
#' to other data tables using either of Sample_Id or Subject_Id.
#'
#' @inheritParams load_wes_mut
#' @return a data.frame with summary
#' @importFrom readr cols
#' @export
load_summary_wes_mut <- function(study, mode="oncotator_filtered"){
  load_from_data(get_filepath_summary_wes_mut(study, mode), progress=F, show_col_types=F)
}
