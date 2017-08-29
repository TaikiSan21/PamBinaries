%%  WARNING: This file was automatically generated from the associated 
%%  sapa_sdf.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the sapa_sdf.mid file in the project MID directory. Once the 
%%  sapa_sdf.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the ACVS function

\name{ACVS}
\alias{ACVS}
\title{Autocovariance sequence}
\concept{utilities}
\usage{ACVS(x, biased=TRUE, center=TRUE)}
\description{Calculates the autocovariance sequence for an input time series.}

\arguments{
\item{x}{a numeric vector representing a uniformly sampled real-valued time series.}
\item{biased}{a logical value. If \code{TRUE}, the biased estimator (normalized by \eqn{N}, the number of samples in the time series) is returned.
If \code{FALSE}, the result is the unbiased estimator (the kth ACVS value is normalized by \eqn{N - |k|} for
the unbiased case where \eqn{k=0,\ldots,N-1}{k=0,...,N-1}). Default: \code{TRUE}.}

\item{center}{a logical value. If \code{TRUE}, the series is first centered (sample mean is subtracted
from series) prior to calculating the ACVS. Default: \code{TRUE}.}
}

\value{
a numeric vector containing the single-sided ACVS for lags \eqn{k=0,\ldots,N-1}{k=0,...,N-1} where \eqn{N} is the length
of the input time series.
}
\seealso{
\code{\link{SDF}}.}

\examples{
## calculate the ACVS for an N(0,1) realization 
plot(seq(0,99), ACVS(rnorm(100)), type="l", lwd=2,
    xlab="lag",ylab="ACVS(rnorm(100))")
ifultools::gridOverlay()
}
\keyword{univar}

