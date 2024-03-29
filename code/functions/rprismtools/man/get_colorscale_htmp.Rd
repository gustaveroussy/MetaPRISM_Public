% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_Plotly.R
\name{get_colorscale_htmp}
\alias{get_colorscale_htmp}
\title{Colorscale vector}
\usage{
get_colorscale_htmp(z, colors_limits = NULL, colors_palette = NULL)

get_colorscale_htmp(z, colors_limits = NULL, colors_palette = NULL)
}
\arguments{
\item{z}{a numeric matrix or vector.}

\item{colors_limits}{(optional) A numeric vector, used to associate to each unique value in \code{z} an interval.}

\item{colors_palette}{(optional) A name ("Reds" or "Blues") or a character vector of color names/codes.}
}
\value{
a \code{data.frame}.

a \code{data.frame}.
}
\description{
Associate to each unique value in the numeric matrix \code{z} a color using a set of colors (a palette) and a set of
limits defining intervals with the same color.

Associate to each unique value in the numeric matrix \code{z} a color using a set of colors (a palette) and a set of
limits defining intervals with the same color.
}
\keyword{internal}
