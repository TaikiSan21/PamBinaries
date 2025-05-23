% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readDLModelData.R
\name{readDLModelData}
\alias{readDLModelData}
\title{Read Deep Learning Model Data}
\usage{
readDLModelData(fid, fileInfo, data, debug = FALSE, ...)
}
\arguments{
\item{fid}{binary file identifier}

\item{fileInfo}{structure holding the file header and module header}

\item{data}{a structure containing standard data}

\item{debug}{logical flag to show more info on errors}

\item{\dots}{Arguments passed to other functions}
}
\value{
a structure containing data from a single object, and a logical
  flag if an error has occurred
}
\description{
Reads binary data stored by the Deep Learning Model module
}
\author{
Taiki Sakai \email{taiki.sakai@noaa.gov}
}
