# @created: 08 Apr 21
# @modified: 27 Apr 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_dsg

#' Load the design table
#'
#' Load the design table for the specified study. Only "prism" is supported.
#'
#' @return a data.frame
#' @param study Name of the cohort used in the naming of the files. Choose "met500", "prism" or "tcga".
#'
#' @author Yoann Pradat
#'
#' @export
load_dsg <- function(study){
  load_from_data(get_filepath_dsg_curated(study), progress=F, show_col_types=F)
}
