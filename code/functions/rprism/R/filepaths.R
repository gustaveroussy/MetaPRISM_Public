
#' Read yaml file of filepath into list.
#' @keywords internal
get_filepaths <- function(){
  cwd <- setwd_to_data()
  filepath <- "organisation/DATA_ORGANISATION.yaml"
  FILEPATHS <- yaml::read_yaml(filepath)
  setwd(cwd)

  FILEPATHS
}

#' Filepath to curated sample data files
#' @param study Name of the study. Choose 'met500', 'prism' or 'tcga'.
#' @param mode Choose the file. 'all' for data about all samples or 'in_design' for data about samples in the design
#' only.
#' @keywords internal
get_filepath_bio_curated <- function(study, mode){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["bio"]][[mode]]
}



#' Filepath to curated subject data files
#' @param study Name of the study. Choose 'met500', 'prism' or 'tcga'.
#' @param mode Choose the file. 'all' for data about all subjects or 'in_design' for data about subjects in the design
#' only.
#' @keywords internal
get_filepath_cln_curated <- function(study, mode){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["cln"]][[mode]]
}


#' Filepath to curated design data table
#'
#' @param study Name of the study. Only 'prism' is supported.
#' @keywords internal
get_filepath_dsg_curated <- function(study){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["dsg"]]
}


#' Filepath to curated id files
#'
#' @param study Name of the study. Choose 'met500', 'prism' or 'tcga'.
#' @keywords internal
get_filepath_ids_curated <- function(study){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["ids"]]
}


#' Filepath to a resource file from a database.
#'
#' @param database Name of the database.
#' @param name Name of the resource.
#' @keywords internal
get_filepath_res <- function(database, name){
  FILEPATHS <- get_filepaths()
  FILEPATHS[["resources"]][[database]][[name]]
}


#' Filepath to colors.
#'
#' @keywords internal
get_filepath_colors <- function(){
  "data_overview/colors/colors.xlsx"
}


#' Filepath to raw mutation tables
#'
#' @inheritParams load_wes_mut
#' @keywords internal
get_filepath_wes_mut <- function(study, mode){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["wes_mut"]][[mode]]
}


#' Filepath to rna-seq gene fusion tables
#'
#' @inheritParams load_rna_fus
#' @keywords internal
get_filepath_rna_fus <- function(study, mode){
  FILEPATHS <- get_filepaths()
  return(FILEPATHS[[study]][["rna_fus"]][[mode]])
}


#' Filepath to rna expression tables
#'
#' @inheritParams load_rna_gex
#' @keywords internal
get_filepath_rna_gex <- function(study, level="genes", metric="counts", test_mode=F, other_mode=NULL){
  FILEPATHS <- get_filepaths()
  if (!level %in% FILEPATHS[[study]][["rna_gex"]]){
    if (!is.null(other_mode)){
      mode <- other_mode
    } else {
      if (test_mode){
        mode <- "test"
      } else {
        mode <- "full"
      }
    }
    return(FILEPATHS[[study]][["rna_gex"]][[mode]][[level]][[metric]][["data"]])
  } else {
    return(FILEPATHS[[study]][["rna_gex"]][[level]][[metric]][["data"]])
  }
}


#' Filepath to summary rna-seq gene expression tables
#'
#' @inheritParams load_rna_gex
#' @keywords internal
get_filepath_summary_rna_gex <- function(study, level="genes", metric="counts", test_mode, other_mode=NULL){
  FILEPATHS <- get_filepaths()
  if (!level %in% FILEPATHS[[study]][["rna_gex"]]){
    if (!is.null(other_mode)){
      mode <- other_mode
    } else {
      if (test_mode){
        mode <- "test"
      } else {
        mode <- "full"
      }
    }
    return(FILEPATHS[[study]][["rna_gex"]][[mode]][[level]][[metric]][["summary"]])
  } else {
    return(FILEPATHS[[study]][["rna_gex"]][[level]][[metric]][["summary"]])
  }
}

#' Filepath to summary rna-seq gene fusion tables
#'
#' @inheritParams load_rna_fus
#' @keywords internal
get_filepath_summary_rna_fus <- function(study, mode){
  FILEPATHS <- get_filepaths()
  return(FILEPATHS[[study]][["rna_fus"]][["summary"]][[mode]])
}

#' Filepath to summary mutation tables
#'
#' @inheritParams load_wes_mut
#' @keywords internal
get_filepath_summary_wes_mut <- function(study, mode){
  FILEPATHS <- get_filepaths()
  FILEPATHS[[study]][["wes_mut"]][["summary"]][[mode]]
}

