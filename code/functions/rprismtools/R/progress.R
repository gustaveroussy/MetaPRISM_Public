#' A progress bar function to run within a \code{for} loop
#'
#' This function should be run within a \code{for} loop, it produces a progress bar to indicate how many iterations have
#' passed 
#' @param i A numeric value indicating the current iteration of the containing for loop
#' @param n A numeric value indicating the total number of iterations of the containing for loop
#' @keywords progress bar for
#' @author Joseph Crispell see repo \url{https://github.com/JosephCrispell/basicPlotteR}
#' @export
#' @examples
#' # Set the number of iterations of the for loop
#' n <- 1000
#' 
#' for(i in 1:n){
#'  
#'   # Sleep for a tenth of a second
#'   Sys.sleep(0.01)
#'   
#'   # Update the progress bar
#'   progress(i, n)
#' }
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
