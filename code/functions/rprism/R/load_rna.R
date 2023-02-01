# @created: 29 Oct 20
# @modified: 30 Aug 21
# @author: Yoann Pradat
# 
# Exported functions
#     load_rna_fus
#     load_rna_gex


#' Load fusions called from rna-seq data.
#'
#' This function loads RNA gene fusion tables for the specified cohort or subcohort.
#'
#' @return a tibble data.frame
#'
#' @param study Name of the cohort used in the naming of the files. Choose 'met500', 'prism' or 'tcga'.
#' @param mode The mode determines the caller and/or the filtering level.
#' \itemize{
#' \item "met500": arriba, ericscript, pizzly, starfusion, oncokb_civic.
#' \item "prism": arriba, ericscript, pizzly, starfusion, oncokb_civic.
#' \item "tcga": aggregated, oncokb_civic.
#' }
#' @param identifiers (optional) If not NULL, return RNA data only for the identifiers values. Used only if
#' \code{identifiers_name} is not NULL.
#' @param identifiers_name (optional) Field name of the identifiers.
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @author Yoann Pradat
#'
#' @export
load_rna_fus <- function(study, mode=NULL, identifiers=NULL, identifiers_name=NULL, ...){
  # load maf file
  cat("-reading aggregated RNA table ... ")
  filepath <- get_filepath_rna_fus(study, mode)
  df_rna <- suppressWarnings(load_from_data(filepath, ...))
  cat("ok!\n")

  # subset if any
  df_rna <- subset_data(df_rna, values=identifiers, col_name=identifiers_name)

  df_rna
}

#' Choose which columns should be loaded.
#'
#' Based on a dataframe with summary data that links column names of the raw rna tables to subject and samples
#' identifiers, this function subsets the summary data on the exact set ofcolumns that will be loaded for the parameters
#' specified.
#'
#' @param df_meta a data.frame with summary data, loaded from \code{\link{load_summary_rna_gex}}
#' @inheritParams load_rna_gex
#' @return a data.frame that contains part of the rows of \code{df_meta}
#'
#' @author Yoann Pradat
#' @keywords internal
choose_columns <- function(df_summary, study, metric, identifiers, identifiers_name){
  metric2prefix = list("TPM"    = "abundance",
                       "counts" = "counts",
                       "length" = "length")

  # select on metric
  if (!is.null(metric) & study!="tcga" & study!="tcga_6_samples_gao"){
    prefix <- metric2prefix[[metric]]
    df_subsummary <- df_summary[grepl(paste0("^", prefix), df_summary[["Col_Name"]]),]
  } else {
    df_subsummary <- df_summary
  }

  # select on identifiers
  if (!is.null(identifiers) & !is.null(identifiers_name)){
    df_subsummary <- df_subsummary[df_subsummary[[identifiers_name]] %in% identifiers,]
  }

  df_subsummary
}


#' Load raw gene expression data files.
#'
#' This function loads RNA gene expression tables for the specified cohort or subcohort.
#'
#' @return a tibble data.frame
#'
#' @param study Name of the cohort used in the naming of the files. Choose one of.
#' \itemize{
#' \item "met500": Trim Galore > Kallisto > Tximport. See \url{https://github.com/gustaveroussy/MetaPRISM}
#' \item "prism": Trim Galore > Kallisto > Tximport. See \url{https://github.com/gustaveroussy/MetaPRISM}
#' \item "tcga": From \url{https://stanfordmedicine.app.box.com/s/lu703xuaulfz02vgd2lunxnvt4mfvo3q}
#' \item "tcga_6_samples_prism": PRISM pipeline on 6 TCGA FASTQs.
#' \item "tcga_6_samples_gao": Gao pipeline on 6 TCGA FASTQs. 
#' }
#' @param level Level at which counts are aggregated. You may choose "genes" or "transcripts".
#' @param metric string or vector of strings. If not specified, every metric is loaded.
#' \itemize{
#' \item "TPM" to load "abundance.kallisto"
#' \item "counts" to load "counts.kallisto"
#' \item "length" to load "length.kallisto"
#'  }
#' @param test_mode (optional) Set to True to use submatrices generated with \code{generate_rna_test} from
#' \code{_util_rna} module.
#' @param identifiers (optional) If not NULL, return RNA data only for the identifiers values. Used only if
#' \code{identifiers_name} is not NULL.
#' @param identifiers_name (optional) Field name of the identifiers.
#' @param ... Extra parameters passed to \code{\link{load_from_data}}
#'
#' @importFrom readr cols_only
#' @importFrom dplyr rename
#'
#' @author Yoann Pradat
#'
#' @export
load_rna_gex <- function(study, level=c("genes", "transcripts"), metric=c("counts", "TPM", "length"), test_mode=F,
                     identifiers=NULL, identifiers_name=NULL, ...){
  level <- match.arg(level)
  metric <- match.arg(metric)

  # load summary data
  df_summary <- load_summary_rna_gex(study, level=level, metric=metric, test_mode=test_mode)

  # choose cols
  df_summary <- choose_columns(df_summary, study, metric, identifiers, identifiers_name)
  usecols <- df_summary[["Col_Name"]]
  col_types <- as.list(rep("d", times=length(usecols)))
  names(col_types) <- usecols

  # load maf file
  cat("-reading aggregated RNA table ... ")
  filepath <- get_filepath_rna_gex(study, level, metric, test_mode)
  if (study=="tcga"){
    col_types[["geneID"]] <- "c"
    col_types <- do.call(readr::cols_only, col_types)
    df_rna <- suppressWarnings(load_from_data(filepath, delim="\t", col_types=col_types, ...))
    df_rna <- dplyr::rename(df_rna, ensembl_gene_id=geneID)
  } else {
    col_types[["ensembl_gene_id"]] <- "c"
    col_types <- do.call(readr::cols_only, col_types)
    df_rna <- load_from_data(filepath, col_types=col_types, ...)
  }
  cat("ok!\n")

  df_rna
}

