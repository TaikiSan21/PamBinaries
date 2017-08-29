%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsDistanceMetric function

\name{mutilsDistanceMetric}
\alias{mutilsDistanceMetric}
\title{L-p metric conversion}
\concept{utilities}
\usage{mutilsDistanceMetric(metric)}
\description{Converts a numeric value (p) representing an L-p metric to an enumerated index
appropriate for Insightful's MUTILS C library, used by many functions in the FRACTAL package
involving \code{itCall} calls.}

\arguments{
\item{metric}{an integer (p) representing an L-p Euclidean distance metric.
Supported values are 1, 2 and \code{Inf}.}

}

\value{
an integer representing the equivalent MUTILS enumerated type for the specifed metric.
}
\seealso{
\code{\link{mutilsFilterType}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsSDF}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
mutilsDistanceMetric(Inf)
}
\keyword{IO}

