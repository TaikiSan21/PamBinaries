%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsFilterTypeContinuous function

\name{mutilsFilterTypeContinuous}
\alias{mutilsFilterTypeContinuous}
\title{Converts a continuous wavelet filter string to MUTILS enum type}
\concept{utilities}
\usage{mutilsFilterTypeContinuous(x)}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function maps a conntinuous wavelet filter character string to the
corresponding enum type in MUTILS.}

\arguments{
\item{x}{a character string representing the continuous wavelet filter type.}

}

\value{
an integer corresponding to the MUTILS enum.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterType}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
## map the Morlet filter 
mutilsFilterTypeContinuous("morlet")
}
\keyword{IO}

