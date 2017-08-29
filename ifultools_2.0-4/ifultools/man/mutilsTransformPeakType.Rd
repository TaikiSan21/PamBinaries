%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsTransformPeakType function

\name{mutilsTransformPeakType}
\alias{mutilsTransformPeakType}
\title{Converts wavelet transform local peak type string to MUTILS enum type}
\concept{utilities}
\usage{mutilsTransformPeakType(x)}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function maps a wavelet transform local peak type character string to the
corresponding enum type in MUTILS.}

\arguments{
\item{x}{a character string representing the wavelet transform local peak type.
Choices are \code{"extrema"}, \code{"maxima"}, and \code{"minima"}.}

}

\value{
an integer corresponding to the MUTILS enum.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterType}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
## map the peak type 
mutilsTransformPeakType("maxima")
}
\keyword{IO}

