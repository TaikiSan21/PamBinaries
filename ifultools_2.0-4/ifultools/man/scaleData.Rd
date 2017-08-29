%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the scaleData function

\name{scaleData}
\alias{scaleData}
\title{Scale numeric data}
\concept{utilities}
\usage{scaleData(x, scale.="linear")}
\description{Scales the input data via a selected logarithmic function.}

\arguments{
\item{x}{an object that inherits from the numeric class, typically a \code{scalar}, \code{vector}, or \code{matrix}.}
\item{scale.}{a character string defining the type of scaling to perform.
Choices are \code{"linear"},\code{"log2"},\code{"log10"},\code{"log"} or \code{"db"}. Default: \code{"linear"} (no scaling).}
}

\value{
the scaled data is returned with an attribute \code{"scalestr"} attached, which defines
the scaling treatment used on the input data.
}
\seealso{
\code{\link[base:exp]{logb}}, \code{\link{decibel}}.}

\examples{
scaleData(c(1,10,100,1000),scale="db")
}
\keyword{utilities}

