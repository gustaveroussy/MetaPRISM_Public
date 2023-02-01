# modified: 19 May 2021
# @created: 06 Nov 2020
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
# 
#     Institut Gustave Roussy
#     Prism Center
#     114 rue Edouard Vaillant, Villejuif, 94800 France
# 
# Fonctions for setting the current working directory to folders of the project_
# 
# The following functions are implemented
#     setwd_to_data()
#     setwd_to_results()
#     setwd_to_logs()
#     setwd_to_scripts()

# Name of the project root folder_ Modify if you've cloned the git repository under a different name
PROJECT_FOLDER <- "MetaPRISM"

#' Set the working directory to the data folder of the project.
#'
#' @return current_wd name of the working directory before change
#' @author Yoann Pradat
#' @export
setwd_to_data <- function(){
  current_wd <- getwd()
  if (!grepl(PROJECT_FOLDER, getwd())){
    stop(paste("Please set the working directory to a location in the repository", PROJECT_FOLDER))
  } else {
    while (!endsWith(getwd(), PROJECT_FOLDER)){
      setwd("..")
    }
  }
  setwd("./data")
  return(current_wd)
}

#' Set the working directory to the results folder of the project.
#'
#' @return current_wd name of the working directory before change
#' @author Yoann Pradat
#' @export
setwd_to_results <- function(){
  current_wd <- getwd()
  if (!grepl(PROJECT_FOLDER, getwd())){
    stop(paste("Please set the working directory to a location in the repository", PROJECT_FOLDER))
  } else {
    while (!endsWith(getwd(), PROJECT_FOLDER)){
      setwd("..")
    }
  }
  setwd("./results")
  return(current_wd)
}

#' Set the working directory to the logs folder of the project.
#'
#' @return current_wd name of the working directory before change
#' @author Yoann Pradat
#' @export
setwd_to_logs <- function(){
  current_wd <- getwd()
  if (!grepl(PROJECT_FOLDER, getwd())){
    stop(paste("Please set the working directory to a location in the repository", PROJECT_FOLDER))
  } else {
    while (!endsWith(getwd(), PROJECT_FOLDER)){
      setwd("..")
    }
  }
  setwd("./logs")
  return(current_wd)
}

#' Set the working directory to the scripts folder of the project.
#'
#' @return current_wd name of the working directory before change
#' @author Yoann Pradat
#' @export
setwd_to_scripts <- function(){
  current_wd <- getwd()
  if (!grepl(PROJECT_FOLDER, getwd())){
    stop(paste("Please set the working directory to a location in the repository", PROJECT_FOLDER))
  } else {
    while (!endsWith(getwd(), PROJECT_FOLDER)){
      setwd("..")
    }
  }
  setwd("./scripts")
  return(current_wd)
}
