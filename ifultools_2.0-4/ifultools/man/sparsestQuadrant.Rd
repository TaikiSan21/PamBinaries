%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the sparsestQuadrant function

\name{sparsestQuadrant}
\alias{sparsestQuadrant}
\title{Estimates sparsest quadrant in data set}
\concept{utilities}
\usage{sparsestQuadrant(x, y=NULL, nquadrant=5)}
\description{Finds the sparsest are in the current plot for placing an object.}

\arguments{
\item{x}{if \code{y} is not \code{NULL}, this argument is a numeric
vector containing the abscissa (x-axis) values for
the current plot. If \code{y=NULL}, this argument is taken as a numeric
vector of ordinate (y-axis) values.}
\item{nquadrant}{an integer denoting the number of quadrants to partition
the space spanned by x and y. For example, if \code{nquadrants=3}, the
x-y space is partitioned into a 3x3 grid of equisized blocks. Default: 5.}

\item{y}{a numeric vector containing the ordinate (y-axis) values
in the current plot. If \code{NULL}, the \code{x} input argument is taken
as the y-axis values while the x-axis values are extracted from
\code{x} via the \code{positions} function. If there are
no positions available, then the sequence \code{1:length(x)} is used
as the x-axis positions. Default: \code{NULL}.}
}

\value{
a list containing the \code{x} and \code{y}-coordinates of the sparsest
data location in the current plot. In addition,
a two-element \code{corner} vector is returned and should be used directly
in placing a key at the returned coordinates via the \code{autoKey} function.
}
\details{
%
Partitions the space spanned by the input x-y data coordinates
into a uniform \code{nquadrant} by \code{nquadrant} grid and finds
the "sparsest quadrant": defined as the quadrant
that contains a minimum population (pmin) of x-y data points
and is and farthest (in an L-inf sense) from quadrants
containing a population greater than pmin.
}
\seealso{
\code{\link{autoKey}}, \code{\link{autoText}}.}

\examples{
## find the sparset quadrant in a circle 
tt <- seq(0,2*pi,length=100)
x  <- sin(tt)
y  <- cos(tt)
sparsestQuadrant(x,y)
}
\keyword{hplot}

