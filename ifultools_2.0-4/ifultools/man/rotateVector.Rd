%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the rotateVector function

\name{rotateVector}
\alias{rotateVector}
\title{Circularly vector rotation}
\concept{utilities}
\usage{rotateVector(x, shift=1)}
\description{Circularly rotates elememnts in a vector.}

\arguments{
\item{x}{a vector.}
\item{shift}{an integer denoting the shift (number of elements to circularly rotate).
A negative value implies a negative (leftward) shift of the data.}
}

\value{
a numeric vector containing the result.
}
\seealso{
\code{\link{isVectorAtomic}}.}

\examples{
rotateVector(1:5, 2)
rotateVector(1:5, -2)
}
\keyword{utilities}

