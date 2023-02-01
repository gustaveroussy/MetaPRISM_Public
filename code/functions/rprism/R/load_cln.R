# @created: 25 Jan 21
# @modified: 27 Apr 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_cln

#' Load clinical files.
#'
#' This function loads clinical tables (one line is one subject id) for the specified study and optionally specified
#' identifiers. The \code{mode} parameters allow to choose which table should be loaded between 'in_design' and 'all'.
#'
#' @return a data.frame
#'
#' @param study Name of the cohort used in the naming of the files. Choose "met500", "prism" or "tcga".
#' @param identifiers (optional) If not NULL, return clinical tables only for the individuals specified.
#' @param identifiers_name (optional) The name of the column the \code{identifiers} values correspond to.
#' @param mode (optional) Either 'in_design' or 'all' for 'pancancer' annotations or 'brca' for per-cancer annotations.
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @author Yoann Pradat
#'
#' @export
load_cln <- function(study, identifiers=NULL, identifiers_name=NULL, mode="in_design", ...){
  df_cln <- load_from_data(get_filepath_cln_curated(study, mode), progress=F, show_col_types=F, ...)
  df_cln <- subset_data(df=df_cln, values=identifiers, col_name=identifiers_name)

  return(df_cln)
}
