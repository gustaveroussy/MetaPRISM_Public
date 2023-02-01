#' Correct id of normal sample.
#'
#' @param x sample id
#' @importFrom stringr str_extract
#' @keywords internal
get_sample_id_from_germline <- function(x){
  if (grepl("_on", x)){
    x <- str_extract(x, "^[\\w\\-]+(?=_on)")
  }
  gsub("^Sample_", "", x)
}


#' Correct id of tumor sample.
#'
#' @param x sample id
#' @keywords internal
get_sample_id_from_somatic <- function(x){
  x <- gsub("^Sample_|_run1$|_run2$", "", x)
  if (x=="MR254_T2-ADN") x <- "MR254-T2-ADN"

  x
}

