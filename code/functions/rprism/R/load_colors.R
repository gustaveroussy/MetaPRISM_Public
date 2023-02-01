
#' Load lists of colors from Excel table
#'
#' Load one or multiple sheets from an excel workbook containing lists of colors for categories of different fields.
#'
#' @return A named list or a \code{tibble} dataframe. Names are categories of some fields or custom names and values are
#' colors in hexadecimal character vectors.
#' @param sheet Name of the excel spreadsheet.
#' @param as_tibble (optional) Set to TRUE if you would like to have the list of colors returned in a \code{tibble}
#' dataframe. If FALSE, colors are returned as a named list.
#' @importFrom readxl read_excel
#'
#' @export
load_colors <- function(sheet, as_tibble=F){
  cwd <- setwd_to_results()
  filepath <- get_filepath_colors()
  df_colors <- read_excel(path=filepath, sheet=sheet)
  setwd(cwd)

  if (as_tibble){
    return(df_colors)
  } else {
    colors <- list()
    for (i in 1:nrow(df_colors)){
      colors[[df_colors[i, "Name", drop=T]]] <- df_colors[i, "Color", drop=T]
    }
    return(colors)
  }
}
