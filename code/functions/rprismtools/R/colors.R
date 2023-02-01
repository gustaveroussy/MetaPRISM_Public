
#' Associate colors to labels
#'
#' Takes as input a vector of labels (may be non unique) and returns a vector of the same size with one color for each
#' unique label in the labels. The set of unique labels may be specified. It is useful if not all labels are present in
#' the vector of labels.
#'
#' @param labels vector.
#' @param pal (optional) name of a palette of RColorBrewer.
#' @param labels_unique (optional) If not NULL, used to defined colors for labels that are not in the labels vector.
#'   Useful when you wish to harmonize colors between plots that do not have all the labels each.
#' @param alpha (optional) double in \[0,1\].
#' @return a character vector of colors of same size as labels.
#'
#' @author Yoann Pradat
#' @export
get_label_colors <- function(labels, pal="Dark2", labels_unique=NULL, alpha=1){
  if (is.null(labels_unique)){
    labels_unique <- sort(unique(labels))
  }
  palette <- RColorBrewer::brewer.pal(n=RColorBrewer::brewer.pal.info[pal, "maxcolors"], pal)

  lab2col <- list()
  i <- 1
  for (lab in labels_unique){
    lab2col[[lab]] <- grDevices::adjustcolor(palette[(i-1) %% length(palette) + 1], alpha)
    i <- i + 1
  }

  colors <- c()
  for (lab in labels){
    colors <- c(colors, lab2col[[lab]])
  }

  colors
}


#' Display a list of colors
#'
#' From a named list of colors, draw a set of rectangles showing the color and the name inside the rectangle.
#'
#' @param line How many lines in the facetted plot.
#' @param col How many columns in the facetted plot.
#' @param colors Named list of colors.
#' @param cex Font size.
#'
#' @importFrom graphics par rect text
#'
#' @author Yoann Pradat
#' @export
rect_plot_colors <- function(line=NULL, col=NULL, colors, cex=1){
  n <- length(colors)

  if (is.null(line) & is.null(col)){
    stop("One of 'col' or 'line' must be not NULL")
  }
  if (is.null(col)){
    col <- ceiling(n/line)
  }
  if (is.null(line)){
    line <- ceiling(n/col)
  }

  par(mar=c(0,0,0,0))

  # Empty chart
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")

  # Add color background
  rect(  
    rep((0:(col - 1)/col),line)[1:n],  
    sort(rep((0:(line - 1)/line),col),decreasing=T)[1:n],   
    rep((1:col/col),line)[1:n], 
    sort(rep((1:line/line),col),decreasing=T)[1:n],  
    border="white",
    col=unlist(colors))

  # Color names
  text(  
    rep((0:(col - 1)/col),line)[1:n] + 1/(2*col),  
    sort(rep((0:(line - 1)/line),col),decreasing=T)[1:n]+1/(2*line), 
    names(colors), 
    cex=cex, adj=0.5)
}
