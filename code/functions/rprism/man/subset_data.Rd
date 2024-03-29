% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_load.R
\name{subset_data}
\alias{subset_data}
\title{Load a table located in the data folder}
\usage{
subset_data(df, values = NULL, col_name = NULL)
}
\arguments{
\item{df}{a data.frame that contains the column \code{col_name}.}

\item{values}{a list of values to be selected.}

\item{col_name}{Name of the the column to be used to susbet the dataframe.}
}
\value{
a data.frame containing a subset of the rows of \code{df}.
}
\description{
Load a table located in the data folder
}
