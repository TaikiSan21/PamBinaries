%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the autoText function

\name{autoText}
\alias{autoText}
\title{Automatic text placement}
\concept{utilities}
\usage{autoText(x, y=NULL, text.="", cex=1, col=1, nquadrant=5)}
\description{Automatically places a given character string in the sparsest
location of the space spanned by the
the x and y input coordinates.}

\arguments{
\item{x}{if \code{y} is not \code{NULL}, this argument is a numeric
vector containing the abscissa (x-axis) values for
the current plot. If \code{y=NULL}, this argument is taken as a numeric
vector of ordinate (y-axis) values.}
\item{cex}{\code{par} function character expansion value. Default: 1.}

\item{col}{\code{par} function color value for the text. Default: 1.}

\item{nquadrant}{an integer denoting the number of quadrants to partition
the space spanned by x and y. For example, if \code{nquadrants=3}, the
x-y space is partitioned into a 3x3 grid of equisized blocks. Default: 5.}

\item{text.}{a character string to place in the current plot.
Default: \code{""}}

\item{y}{a numeric vector containing the ordinate (y-axis) values
in the current plot. If \code{NULL}, the \code{x} input argument is taken
as the y-axis values while the x-axis values are extracted from
\code{x} via the \code{positions} function. If there are
no positions available, then the sequence \code{1:length(x)} is used
as the x-axis positions. Default: \code{NULL}.}
}

\value{
no output is returned.
}
\note{The defined text is positioned in the estimated sparsest
region of the specified data set.}

\seealso{
\code{\link[graphics]{text}}, \code{\link{autoKey}}, \code{\link{sparsestQuadrant}}, \code{\link{mergeList}}.}

\examples{
## quadratic 
x <- seq(-5,5)
y <- x^2
plot(x,y,type="p")
autoText(x, y, text="Some text in a parabola", cex=1.5)
}
\keyword{hplot}

