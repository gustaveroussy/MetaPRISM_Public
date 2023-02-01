# @created: 29 Oct 20
# @modified: 19 Feb 21
# @author: Yoann Pradat
#
# Internal functions
#     progress
# Exported functions
#     merge_on_rows

# Internal =============================================================================================================

#' A progress bar function to run within a \code{for} loop
#'
#' This function should be run within a \code{for} loop, it produces a progress bar to indicate how many iterations have
#' passed 
#' @param i A numeric value indicating the current iteration of the containing for loop
#' @param n A numeric value indicating the total number of iterations of the containing for loop
#' @author Joseph Crispell see repo \url{https://github.com/JosephCrispell/basicPlotteR}
#'
#' @keywords internal
progress <- function(i, n){
  # Note the width of the console
  consoleWidth <- options("width")$width
  progressBarWidth <- consoleWidth - 10 # 10 characters to leave space for dial, bar ends, and percentage
  
  # Calculate the percentage done
  percentage <- round((i / n)*100, digits=0)
  if(nchar(percentage) != 3){
    percentage <- paste0(paste(rep(" ", 3-nchar(percentage)), collapse=""), percentage)
  }
  
  # Check if on first iteration - create the progress bar
  if(i == 1){
    cat(paste0("|", paste(rep(" ", progressBarWidth), collapse=""), "| | ", percentage, "%"))
  
  # Update the progress bar
  }else if(i != n){
    
    # Calculate how many progress bar parts we have progressed
    nParts <- ceiling(i / (n/progressBarWidth))
    
    # Create the tenths bars and pad with spaces for the ones that are absent
    progressString <- paste0(paste(rep("-", nParts), collapse=""), 
                             paste(rep(" ", progressBarWidth - nParts), collapse=""))

    # Add spinning dial
    remainderFromFour <- i %% 4
    if(remainderFromFour == 0){
      dial <- "|"
    }else if(remainderFromFour == 1){
      dial <- "/"
    }else if(remainderFromFour == 2){
      dial <- "\u2500"
    }else if(remainderFromFour == 3){
      dial <- "\\"
    }

    # Update the progress bar
    cat(paste0("\r|", progressString, "| ", dial, " ", percentage, "%"))
  
  # Add a finished statement
  }else{
    cat(paste0("\r|", paste0(rep("-", progressBarWidth), collapse=""), "|   100%\n"))
  }
}

# Export ===============================================================================================================

#' Merge two dataframes on their row names.
#'
#' This function merges two dataframes on their row names and provides the user with the choice on which columns should 
#' be kept.
#'
#' @return a data.frame
#'
#' @param df_x a data.frame
#' @param df_y a data.frame
#' @param how_cols a character vector specifying how the merge should proceed
#' \enumerate{
#' \item 'inner' keep only the columns common to both \code{df_x} and \code{df_y}
#' \item 'x' keep only the columns of \code{df_x} 
#' \item 'y' keep only the columns of \code{df_y} 
#' \item 'outer' keep all the columns of \code{df_x} and \code{df_y} 
#' }
#'
#' This function raises an error if a row entry (identified by its row name) is present in both \code{df_x} and
#' \code{df_y} and has divergent values on a column that is present in both data frames.
#'
#' @author Yoann Pradat
#'
#' @keywords export
merge_on_rows <- function(df_x, df_y, how_cols="inner"){
  how_cols <- match.arg(how_cols, choices=c("inner", "x", "y", "outer"))

  cols_x <- colnames(df_x)
  cols_y <- colnames(df_y)
  cols_i <- intersect(cols_x, cols_y)

  df_m <- merge(df_x, df_y, by=c("row.names", cols_i), all=T)

  if (nrow(df_m) > length(unique(c(rownames(df_x), rownames(df_y))))){
    stop("There is a entry that has the same row_name in df_x and df_y but divergent values")
  }

  rownames(df_m) <- df_m[["Row.names"]]
  df_m["Row.names"] <- NULL

  if (how_cols == "inner"){
    return(df_m[,cols_i])
  } else if (how_cols == "x"){
    return(df_m[,cols_x])
  } else if (how_cols == "y"){
    return(df_m[,cols_y])
  } else {
    return(df_m)
  }
}
