# @created: 21 May 21
# @modified: 21 May 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_resource

#' Load resources from different databases or used for different analyses.
#'
#' @param database Name of the database to be load from. Supported values are
#' \itemize{
#'   \item{civic}
#'   \item{oncokb}
#'   \item{cosmic}
#'   \item{gencode}
#'   \item{curated}
#' }
#' @param name Name of the resource
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @author Yoann Pradat
#' @export
load_resource <- function(database, name, ...){
  df <- load_from_data(get_filepath_res(tolower(database), tolower(name)), ...)
  df
}
