# @created: 30 Dec 21
# @modified: 30 Dec 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_ids

#' Load identifiers files.
#'
#' This function loads the table linking all subject, sample and idspsy identifiers..
#'
#' @return a data.frame
#'
#' @param study Name of the cohort used in the naming of the files. Choose "met500", "prism" or "tcga".
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @author Yoann Pradat
#'
#' @export
load_ids <- function(study, ...){
  df_ids <- load_from_data(get_filepath_ids_curated(study), progress=F, show_col_types=F, ...)

  df_ids
}
