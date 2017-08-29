%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the stringSize function

\name{stringSize}
\alias{stringSize}
\title{Size of a character string in current font}
\concept{utilities}
\usage{stringSize(x, srt=0, adj=0.5, cex=1)}
\description{Calculates the approximate relative x and
y spans of the input character string x
in the current par("usr") coordinates. The \code{srt},
\code{adj}, and \code{cex} inputs are exactly those specified
by the \code{par} function.}

\arguments{
\item{x}{a character string.}
\item{adj}{a numeric value on the interval [0,1] used to justify the text relative to
its placement in a plot. Default: 0.5 (center aligned).}

\item{cex}{the expansion value for the character string. Default: 1.}

\item{srt}{a numeric value defining the amount of string rotation in degrees. Default: 0.}
}

\value{
a list of x- and y-vectors that can be used to estimate the
width and height, respectively, of the input character string in the current
par("usr") coordinates. Each vector contains a relative starting and ending
value and the absolute difference between the two values repesents the span
in the corresponding direction.
}
\note{A graphics window will open up if one is not already opened.}

\seealso{
\code{\link[graphics]{par}}, \code{\link{autoKey}}, \code{\link{autoText}}, \code{\link{em}}.}

\examples{
stringSize("What's for dinner?", adj=0, srt=45, cex=1.5)
}
\keyword{hplot}

