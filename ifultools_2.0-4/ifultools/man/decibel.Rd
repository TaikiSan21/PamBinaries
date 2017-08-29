%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the decibel function

\name{decibel}
\alias{decibel}
\title{Conversion to decibels}
\concept{utilities}
\usage{decibel(x, type=1, na.zero=TRUE)}
\description{Convert numeric input into decibels.}

\arguments{
\item{x}{a real positive numeric value.}
\item{na.zero}{logical value controlling the handing of zeros in the data. If \code{TRUE}, the
zero values are replaced with NAs prior to converting to decibels, avoiding an error when
attempting to take the logarithm of zero.}

\item{type}{an integer signifying the decible form to use in the conversion.
If \code{type=1} \eqn{10 * log10(x)} is used to convert the input
and \eqn{20 * log10(x)} isued for \code{type=2}.}
}

\value{
the decibel equivalent of the input \code{x}.
}
\seealso{
\code{\link{nextDyadic}}, \code{\link{isVectorAtomic}}.}

\examples{
decibel(c(1,10,100,1000), type=1)
decibel(c(1,10,100,1000), type=2)
}
\keyword{utilities}

