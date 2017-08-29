%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the logScale function

\name{logScale}
\alias{logScale}
\title{Logarithmically spaced scale vector generation}
\concept{utilities}
\usage{logScale(scale.min, scale.max,
    scale.ratio=2, scale.res=NULL, coerce=NULL)}
\description{Creates a vector whose values are base-2 logarithmically spaced.}

\arguments{
\item{scale.min}{an integer denoting the minimum scale.}

\item{scale.max}{an integer denoting the maximum scale.}
\item{coerce}{the coercion function to use on the resulting scale vector. Default: NULL (no coercion).}

\item{scale.ratio}{a numeric value representing the ratio of the next scale over the current scale, i.e., the ratio of successive scales.
Default: 2.}

\item{scale.res}{the numeric value denoting scale resolution. This input serves
as an alternative to the \code{scale.ratio} input, and is related via
\code{scale.ratio <- 1 / scale.res + 1} or conversely \code{scale.res <- 1/(scale.ratio - 1)}.
Default: \code{NULL} (\code{scale.ratio} is used instead).}
}

\value{
a numeric vector containing the scales.
}
\seealso{
\code{\link[base]{scale}}, \code{\link{ilogb}}.}

\examples{
logScale(scale.min = 1, scale.max=34, scale.ratio=4)
}
\keyword{utilities}

