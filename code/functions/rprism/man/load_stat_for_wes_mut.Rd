% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_wes.R
\name{load_stat_for_wes_mut}
\alias{load_stat_for_wes_mut}
\title{Build a dataframe with per patient statistics of MAF files loading.}
\usage{
load_stat_for_wes_mut(
  df_maf,
  patient_field,
  tumor_field = NULL,
  normal_field = NULL
)
}
\arguments{
\item{df_maf}{The dataframe in MAF format}

\item{patient_field}{The name of the patient column}

\item{tumor_field}{The name of the tumor column}

\item{normal_field}{The name of the normal column}
}
\value{
a list of \code{data.frame}
}
\description{
Used to see which data was loaded for which patient in MAF files.
}
\author{
Yoann Pradat
}
