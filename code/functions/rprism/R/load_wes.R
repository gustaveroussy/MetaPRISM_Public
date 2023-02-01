# @created: 25 Jan 21
# @modified: 30 Aug 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_stat_for_wes_mut
#     load_wes_mut

#' Build a dataframe with per patient statistics of MAF files loading.
#'
#' Used to see which data was loaded for which patient in MAF files.
#'
#' @return a list of \code{data.frame}
#'
#' @param df_maf The dataframe in MAF format
#' @param patient_field The name of the patient column
#' @param tumor_field The name of the tumor column
#' @param normal_field The name of the normal column
#'
#' @importFrom dplyr group_by summarize n
#' @importFrom rlang := .data
#' @importFrom magrittr %>%
#'
#' @author Yoann Pradat
#'
#' @export
load_stat_for_wes_mut <- function(df_maf, patient_field, tumor_field=NULL, normal_field=NULL){
  if (is.null(tumor_field) & !is.null(normal_field)){
    df_stat <- unique(df_maf[,c(patient_field, normal_field)])
    rownames(df_stat) <- seq(1, nrow(df_stat))

    df_stat <- df_stat %>%
      dplyr::group_by(.data[[patient_field]]) %>%
      dplyr::summarize(N_Loaded=dplyr::n(), {{ normal_field }}:=paste(.data[[normal_field]], collapse="|"))

    message <- "normal sample"
  } else {
    df_stat <- unique(df_maf[,c(patient_field, tumor_field, normal_field)])
    rownames(df_stat) <- seq(1, nrow(df_stat))

    tumor_vs_normal <- paste(tumor_field, normal_field, sep="_vs_")
    df_stat[, tumor_vs_normal] <- paste(df_stat[[tumor_field]], df_stat[[normal_field]], sep="_vs_")
    df_stat <- df_stat[, c(patient_field, tumor_vs_normal)]

    df_stat <- df_stat %>%
      dplyr::group_by(.data[[patient_field]]) %>%
      dplyr::summarize(N_Loaded=dplyr::n(), {{ tumor_vs_normal }}:=paste(.data[[tumor_vs_normal]], collapse="|"))

    message <- "couple(s) tumor vs normal"
  }

  # print
  s_n_loaded <- table(df_stat$N_Loaded)

  for (n in sort(names(s_n_loaded))){
    if (as.integer(n) <= 2){
      cat(paste(n, message, "loaded for", s_n_loaded[[n]], "patient(s)\n"))
    } else {
      cat(paste(n,  message, "loaded for", s_n_loaded[[n]], "patient(s)\n"))
      cat(paste0("\t", paste(df_stat[df_stat$N_Loaded == as.integer(n), patient_field, drop=T], collapse="\n\t"), "\n"))
    }
  }

  df_stat
}


#' Load mutation tables.
#'
#' Load mutation tables for the specified study and optionally specified identifiers. This function is an R
#' reimplementation of the function \code{load_wes_mut} from the pyprism package.
#'
#' \itemize{
#'  \item If "somatic_maf", load data/[study]/wes/somatic_maf/somatic_calls.maf.gz
#'  \item If "somatic_filters", load data/[study]/wes/somatic_maf/somatic_calls_filters.tsv.gz
#'  \item If "somatic_oncokb", load data/[study]/wes/somatic_maf/somatic_calls_oncokb.maf.gz
#'  \item If "somatic_civic", load data/[study]/wes/somatic_maf/somatic_calls_civic.maf.gz
#' }
#'
#' @return a \code{data.frame}
#'
#' @param study Name of the cohort used in the naming of the files. Choose one of "met500", "prism" or "tcga".
#' @param identifiers (optional) If not NULL, return clinical tables only for the individuals specified
#' @param identifiers_name (optional) The name of the column the \code{identifiers} values correspond to.
#' @param mode (optional) Choose which MAF file to load.
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @importFrom readr cols
#' @importFrom dplyr left_join
#'
#' @author Yoann Pradat
#'
#' @export
load_wes_mut <- function(study, identifiers=NULL, identifiers_name=NULL, mode="somatic_maf", ...){
  # load maf file
  cat("-reading aggregated MAF file ... ")
  df_maf <- load_from_data(get_filepath_wes_mut(study, mode), header_prefix="##", col_types=readr::cols(), ...)
  cat("ok!\n")

  # use summary data to add sample and subject identifiers
  # df_summary <- load_summary_wes_mut(study, mode)
  # df_summary[, "N_Mutations"] <- NULL
  # cols_com <- intersect(colnames(df_summary), colnames(df_maf))
  # df_maf <- left_join(df_maf, df_summary, by=cols_com)

  # subset if any
  df_maf <- subset_data(df_maf, values=identifiers, col_name=identifiers_name)

  df_maf
}
