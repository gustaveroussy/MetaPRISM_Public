% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_summary.R
\name{load_summary_wes_mut}
\alias{load_summary_wes_mut}
\title{Load summary of MAF tables}
\usage{
load_summary_wes_mut(study, mode = "oncotator_filtered")
}
\arguments{
\item{study}{Name of the cohort used in the naming of the files. Choose one of "met500", "prism" or "tcga".}

\item{mode}{(optional) Choose which MAF file to load.}
}
\value{
a data.frame with summary
}
\description{
The summary table links original names in the WES table to Sample_Id and Subject_Id. This allows easier connection
to other data tables using either of Sample_Id or Subject_Id.
}
