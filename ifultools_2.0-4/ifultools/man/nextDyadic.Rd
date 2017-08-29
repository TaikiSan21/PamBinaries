%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the nextDyadic function

\name{nextDyadic}
\alias{nextDyadic}
\title{Next integer power of 2}
\concept{utilities}
\usage{nextDyadic(x)}
\description{Find smallest power of two that is greater than or equal to an input x.}

\arguments{
\item{x}{a real number. May be positive, zero, or negative.}

}

\value{
returns the next power of two that is greater than, or equal to,
the input x. If x is negative, 1 is returned.
}
\seealso{
\code{\link{isVectorAtomic}}, \code{\link{decibel}}.}

\examples{
nextDyadic(15)
}
\keyword{utilities}

