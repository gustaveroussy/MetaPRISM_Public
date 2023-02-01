# @created: 29 Oct 20
# @modified: Apr 27 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_bio

#' Load biospecimen files.
#'
#' This function loads biospecimen tables (one line is one sample id) for the specified study and optionally specified
#' identifiers. The \code{mode} parameters allow to choose which table should be loaded between 'in_design' and 'all'.
#'
#' @return a data.frame
#'
#' @param study Name of the cohort used in the naming of the files. Choose "met500", "prism" or "tcga".
#' @param identifiers (optional) If not NULL, return clinical tables only for the individuals specified
#' @param identifiers_name (optional) The name of the column the \code{identifiers} values correspond to.
#' @param mode (optional) Choose "in_design" to load only data for samples in the design or 'all' to load data for all
#' samples of the cohort.
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @author Yoann Pradat
#'
#' @export
load_bio <- function(study, identifiers=NULL, identifiers_name=NULL, mode="in_design", ...){
  df_bio <- load_from_data(get_filepath_bio_curated(study, mode), progress=F, show_col_types=F, ...)
  df_bio <- subset_data(df=df_bio, values=identifiers,  col_name=identifiers_name)

  df_bio
}
