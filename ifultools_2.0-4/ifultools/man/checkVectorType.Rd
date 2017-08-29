%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the checkVectorType function

\name{checkVectorType}
\alias{checkVectorType}
\title{Check for vector type and class}
\concept{utilities}
\usage{checkVectorType(x, isType="numeric")}
\description{Tests the input for being a vector (as defined
by the \code{isVectorAtomic} function and checks the
class of the vector. This function is meant to be used
internally for function input argument verification
and typically will not be used otherwise.}

\arguments{
\item{x}{an S-PLUS object.}
\item{isType}{a character string defining the class of the
object to be checked ala \code{is(x,isType)}. Default: \code{"numeric"}.}
}

\value{
no output is returned. If these tests fail, a \code{stop} condition
is executed.
}
\seealso{
\code{\link{isVectorAtomic}}, \code{\link{checkScalarType}}, \code{\link{checkRange}}.}

\examples{
## expect no output from the following calls 
checkVectorType(rnorm(1024),"numeric")
checkVectorType(1:3,"integer")
checkVectorType(letters,"character")
}
\keyword{utilities}

