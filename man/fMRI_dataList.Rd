% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dfConn.R
\docType{data}
\name{fMRI_dataList}
\alias{fMRI_dataList}
\title{Sample fMRI time series data matrices}
\format{A list consisting of two matrices, each with 278 rows and 750 columns:
\describe{
  \item{subject1}{sample time series data matrix for subject 1.}
  \item{subject2}{sample time series data matrix for subject 2.}
  
  
  For each data matrix in the list:
  \item{row}{index of regions of interest .}
  \item{column}{fMRI time series data at each time point.}
}}
\usage{
data(fMRI_dataList)
}
\description{
Sample functional magnetic resonance imaging (fMRI) data contains a list of data matrices with consisting time series data of two subjects, with 6 scans for each subject in total.
}
\keyword{datasets}
