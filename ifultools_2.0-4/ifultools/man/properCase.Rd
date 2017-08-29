%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the properCase function

\name{properCase}
\alias{properCase}
\title{Case conversion ala a proper noun}
\concept{utility}
\usage{properCase(x, pre="", post="", sep="")}
\description{Converts a character string with arbitrary
case to a character string whose first character
is capitalized and remaining characters are set to
lower case.}

\arguments{
\item{x}{a character string or vector of character strings.}
\item{post}{a character string to append to each character string
in \code{x} after conversion. Default: \code{""}.}

\item{pre}{a character string to prepend to each character string
in \code{x} after conversion. Default: \code{""}.}

\item{sep}{a character string that is pasted between
the other portions of the combined string ala
pre-sep-formatted x-sep-post. Default: \code{""}.}
}

\value{
a vector of character strings containing the result. The length
of this vector is equal to \code{length(x)}.
}
\seealso{
\code{\link{ordinal}}.}

\examples{
## produces "myBasketBall" "mySoccerBall" 
## "myFootBall" 
properCase(c("baSkeT","SoCcEr","foot"),post="Ball",pre="my")
}
\keyword{utilities}

