%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsFilterType function

\name{mutilsFilterType}
\alias{mutilsFilterType}
\title{Converts wavelet filter string to MUTILS enum type}
\concept{utilities}
\usage{mutilsFilterType(wavelet="s8")}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function maps a wavelet filter character string to the
corresponding enum type in MUTILS.}

\arguments{

\item{wavelet}{a character string representing the wavelet filter type. Default: \code{"s8"}.}
}

\value{
an integer corresponding to the MUTILS enum.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
## map the Haar filter 
mutilsFilterType("haar")
}
\keyword{IO}

