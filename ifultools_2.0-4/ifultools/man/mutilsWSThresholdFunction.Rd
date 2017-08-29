%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsWSThresholdFunction function

\name{mutilsWSThresholdFunction}
\alias{mutilsWSThresholdFunction}
\title{Converts wavelet shrinkage theshold function to MUTILS enum type}
\concept{utilities}
\usage{mutilsWSThresholdFunction(x)}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function maps a wavelet transform character string to the
corresponding enum type in MUTILS.}

\arguments{
\item{x}{a character string representing the waveshrink function.
Choices are \code{"universal"}, \code{"minimax"}, and \code{"adaptive"}.}

}

\value{
a list containg objects \code{index} and \code{threshfun},
containing an integer corresponding to the MUTILS enum and name
of the waveshrink threshold function, respectively.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterType}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}.}

\examples{
## obtain MUTILS enum map for the adaptive 
## threshold function 
mutilsWSThresholdFunction("adaptive")
}
\keyword{IO}

