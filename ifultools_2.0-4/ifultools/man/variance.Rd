%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the variance function

\name{variance}
\alias{variance}
\title{Variance computation, usable in both S-PLUS and R}
\concept{utilities}
\usage{variance(x, na.rm=TRUE, unbiased=FALSE)}
\description{A hybrid of S-PLUS and R definitions of the variance function.}

\arguments{
\item{x}{a numeric \code{vector}, \code{matrix}, or \code{data.frame}.}

\item{na.rm}{a logical flag. If \code{TRUE}, NA values are removed prior to computation. Default: \code{TRUE}.}
\item{unbiased}{a logical flag. If \code{TRUE}, the unbiased sample variance is returned. Default: \code{FALSE}.}
}

\value{
the sample variance of the input.
}
\details{
%
Allows the user to calculate the biased/unbiased variance of the input
and provides the option to eliminate NA values as well.
}
\seealso{
\code{\link[stats]{var}}.}

\examples{
set.seed(100)
x <- rnorm(100)
variance(x, unbiased=TRUE)
variance(x, unbiased=FALSE)
variance(c(x,rep(NA,30)), na.rm=TRUE)
}
\keyword{utilities}

