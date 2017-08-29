%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsTransformType function

\name{mutilsTransformType}
\alias{mutilsTransformType}
\title{Converts wavelet transform string to MUTILS enum type}
\concept{utilities}
\usage{mutilsTransformType(x)}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function maps a wavelet transform character string to the
corresponding enum type in MUTILS.}

\arguments{
\item{x}{a character string representing the wavelet transform type.
Choices are \code{"modwt"}, \code{"modwpt"}, \code{"dwt"}, and \code{"dwpt"}.}

}

\value{
an integer corresponding to the MUTILS enum.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterType}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
## map the MODWPT 
mutilsTransformType("modwpt")
}
\keyword{IO}

