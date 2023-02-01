#' Delineate treatments from compact format to more user-friendly format
#'
#' Compressed data must consist of a name with additional data enclosed in brackets and separated by a common delimiter. 
#'
#' @param df a data.frame that must contain the column \code{col_name}.
#' @param col_name Name of the column to be delineated.
#' @param format Character vector of delineated column names corresponding the format of the compressed data.
#' @param separate_rows (optional) Should data be first unnested on the "|" separator"?
#' @param prefix (optional) Prefix to the new columns that contain the delineated data.
#' @param suffix (optional) Suffix to the new columns that contain the delineated data.
#' @param delineate_dates (optional) Used only if "Dates" is in \code{format}. Split "Dates" into "Date_Min" and
#'   "Date_Max".
#' @param convert_dates_to_datetime (optional) Used only if delineate_dates=True. Should dates be converted to datetime
#'   type?
#' @param sep (optional) The delimiter of the compressed data.
#' @param name (optional) Name given to the column containing the value before the brackets. Prefix and suffixes
#' specified will be added to this name.
#' @return a data.frame with possibly more rows (if \code{separate_rows} is set to \code{TRUE}) and more columns
#'   (as many as names in the list \code{format}) built from \code{df}.
#'
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom tidyr separate_rows separate
#'
#' @export
delineate_data <- function(df, col_name, format, separate_rows=T, prefix="", suffix="", delineate_dates=F,
                           convert_dates_to_datetime=T, sep=";", name="Name"){ 

  if (separate_rows){
    df <- df %>% 
      tidyr::separate_rows({{ col_name }}, sep="\\|")
  }

  col_delineate <- "data"
  regex_data <- "(?<=\\()(.*)(?=\\))"

  extract_name <- function(s){
    regex_name <- "^(.*)(?=\\()"
    regex_out <- regexpr(regex_name, s, perl=T)
    if (is.na(regex_out)){
      return(NA)
    } else {
      return(regmatches(s, regex_out))
    }
  }

  df[[name]] <- sapply(df[[col_name]], extract_name, simplify=T, USE.NAMES=F)
  df[[col_delineate]] <- sapply(df[[col_name]], function(s) regmatches(s, regexpr(regex_data, s, perl=T)), USE.NAMES=F)

  format_names <- paste0(prefix, format, suffix)
  suppressWarnings({df <- df %>%
    tidyr::separate(.data[[col_delineate]], format_names, sep=sep)})
  df[[format_names]] <- df[[format_names]] %>% replace(.=="character(0)", NA)

  if (delineate_dates & "Dates" %in% format){
    col_dates <- paste0(prefix, "Dates", suffix)
    col_date_min <- paste0(prefix, "Date_Min", suffix)
    col_date_max <- paste0(prefix, "Date_Min", suffix)
    df <- df %>%
      tidyr::separate({{ col_dates }}, c(col_date_min, col_date_max), sep="_to_")
  }

  if (convert_dates_to_datetime){
    cols_date <- format_names[grepl("Date", format_names) & !grepl("Dates", format_names)]
    df <- df %>%
      dplyr::mutate_at(cols_date, as.Date)
  }

  df
}


#' Load header of a table.
#'
#' @param path Path to the file.
#' @param prefix Prefix for identifying header
#'
#' @keywords internal
read_header <- function(path, prefix){
  hea <- list()
  if (grepl(".gz$", path)){
    con <- gzfile(path, "rt")
  } else {
    con <- file(path, "r")
  }
  i <- 1
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 | !grepl(paste0("^", prefix), line)) {
      break
    } else {
      hea[[i]] <- line
      i <- i+1
    }
  }
  close(con)

  hea
}


#' Load a table.
#'
#' @param path Path to the file.
#' @param header_prefix If not NULL, skip header rows that will be identified as first lines starting with this prefix.
#' @param ... Extra parameters passed to \code{read_delim}.
#' 
#' @importFrom readr read_delim
#' @importFrom readxl read_excel
#'
#' @export
load_table <- function(path, header_prefix=NULL, ...){
  if (is.null(path)) return(NULL)
  args <- list(...)

  if (!"delim" %in% names(args)){
    if (grepl(".csv", path)){
      args$delim <- ","
    } else if (grepl(".tsv", path) | grepl(".txt", path)){
      args$delim <- "\t"
    }
  }

  if (!is.null(header_prefix)){
    header <- read_header(path, prefix=header_prefix)
    args$skip <- length(header) 
  }

  if (grepl(".xlsx$", path)){
    df <- do.call(readxl::read_excel, c(path, args))
  } else {
    if (!"progress" %in% names(args)){
      args$progress <- F
    }
    if (!"show_col_types" %in% names(args)){
      args$show_col_types <- F
    }

    if(grepl(".gz$", path)) {
      file <- base::gzfile(path)
      df <- do.call(readr::read_delim, c(list(file=file), args))
    } else {
      df <- do.call(readr::read_delim, c(path, args))
    }
  }

  df
}

#' Load a table located in the data folder
#'
#' @param path Path to the file.
#' @param ... Extra parameters passed to \code{read_delim}.
#' 
#' @importFrom readr read_delim
#' @importFrom readxl read_excel
#'
#' @export
load_from_data <- function(path, ...){
  cwd <- setwd_to_data()
  df <- load_table(path, ...)
  setwd(cwd)

  df
}


#' Write a table to a path in the data folder
#'
#' @param df a data.frame that contains the column \code{col_name}.
#' @param path a path relative to the folder \code{data/}.
#' @param col_name Name of the the column to be used to susbet the dataframe.
#' @param ... Extra parameters passed to \code{write_xlsx} or \code{write_delim}.
#'
#' @importFrom readr write_delim 
#' @importFrom writexl write_xlsx
#'
#' @export
save_to_data <- function(df, path, ...){
  cwd <- setwd_to_data()
  if (grepl(".xlsx$", path)){
    writexl::write_xlsx(df, path, ...)
  } else {
    if (grepl(".tsv$", path)){
      readr::write_delim(df, path, delim="\t", ...)
    } else if (grepl(".csv$", path)){
      readr::write_delim(df, path, delim=",", ...)
    } else {
      readr::write_delim(df, path, ...)
    }
  }
  cat(paste("-file saved at", path, "\n"))
  setwd(cwd)
}


#' Load a table located in the data folder
#'
#' @param df a data.frame that contains the column \code{col_name}.
#' @param values a list of values to be selected.
#' @param col_name Name of the the column to be used to susbet the dataframe.
#' @return a data.frame containing a subset of the rows of \code{df}.
#'
#' @export
subset_data <- function(df, values=NULL, col_name=NULL){
  if (is.null(col_name) | is.null(values)){
    return(df)
  } else {
    return(df[df[[col_name]] %in% values, ,drop=F])
  }
}
