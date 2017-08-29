%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the checkRange function

\name{checkRange}
\alias{checkRange}
\title{Check the range of a numeric object}
\concept{utilities}
\usage{checkRange(x, range.=0:1, inclusion=c(TRUE,TRUE))}
\description{Tests the input range based on the specified interval.}

\arguments{
\item{x}{an object belonging to class numeric.}
\item{inclusion}{a two-element logical vector defining the boundary tests.
For example, if \code{range.=0:1}, the interval boundaries are tested as follows
\describe{
\item{inclusion=c(TRUE,TRUE)}{range(x) must be on [0,1].}

\item{inclusion=c(TRUE,FALSE)}{range(x) must be on [0,1).}

\item{inclusion=c(FALSE,TRUE)}{range(x) must be on (0,1].}

\item{inclusion=c(FALSE,FALSE)}{range(x) must be on (0,1).}}

If only a single logical element is specified, it is replicated
to form a two-element vector. If more than two elements are specified,
only the first two are used as described above.
Default: \code{c(TRUE,TRUE)}.}

\item{range.}{a two-element numeric vector containing the limits of the test interval.
Default: \code{0:1}.}
}

\value{
no output is returned. If these tests fail, a \code{stop} condition
is executed.
}
\seealso{
\code{\link{isVectorAtomic}}, \code{\link{checkVectorType}}, \code{\link{checkScalarType}}.}

\examples{
## expect no output from the following calls 
checkRange(pi,range=c(0,5))
checkRange(1:5,range=c(1,5), inclusion=c(TRUE,TRUE))
}
\keyword{utilities}

