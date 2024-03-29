% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_load.R
\name{save_to_data}
\alias{save_to_data}
\title{Write a table to a path in the data folder}
\usage{
save_to_data(df, path, ...)
}
\arguments{
\item{df}{a data.frame that contains the column \code{col_name}.}

\item{path}{a path relative to the folder \code{data/}.}

\item{...}{Extra parameters passed to \code{write_xlsx} or \code{write_delim}.}

\item{col_name}{Name of the the column to be used to susbet the dataframe.}
}
\description{
Write a table to a path in the data folder
}
