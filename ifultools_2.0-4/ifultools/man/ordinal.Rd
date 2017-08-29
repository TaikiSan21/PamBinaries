%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the ordinal function

\name{ordinal}
\alias{ordinal}
\title{Ordinal value conversion}
\concept{utilities}
\usage{ordinal(x)}
\description{Converts a given number to the nearest ordinal value.}

\arguments{
\item{x}{a numeric scalar.}

}

\value{
a character string denoting the nearest ordinal value.
}
\seealso{
\code{\link{properCase}}.}

\examples{
ordinal(1)
ordinal(2)
ordinal(3)
ordinal(pi)
ordinal(123)
ordinal(124)
ordinal(1000)
}
\keyword{utilities}

