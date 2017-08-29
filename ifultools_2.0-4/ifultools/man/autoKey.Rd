%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the autoKey function

\name{autoKey}
\alias{autoKey}
\title{Automatic key placement}
\concept{utilities}
\usage{autoKey(x, y=NULL, args.=NULL, nquadrant=5)}
\description{Automatically places a key in the sparsest
location of the space spanned by the
the x and y input coordinates.}

\arguments{
\item{x}{if \code{y} is not \code{NULL}, this argument is a numeric
vector containing the abscissa (x-axis) values for
the current plot. If \code{y=NULL}, this argument is taken as a numeric
vector of ordinate (y-axis) values.}
\item{args.}{a list of arguments that define the legend ala the \code{key} function.
Default: \code{NULL} (no matpoints).}

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
no output is returned.
}
\note{The key as described by \code{args.} is positioned in estimated sparsest
region of the specified data set.}

\details{
%
The argument \code{args.} is a list
of arguments that define the key and \code{nquadrant}
is the number of equi-sized quadrants to use in dividing
the space in both the x- and y-direction. e.g., if
nquadrant=3, the x-y plane is partitioned into a
3x3 uniform grid and the position of the sparsest
quadrant (as defined by \code{sparsetQuadrant})
is used to place the key. This function
does not take into consideration lines connecting
the data points and only considers the spatial distribution
of the points in the x-y plane. Futhermore, this function
assumes that \code{par("usr")} is approximately \code{c(range(x),range(y))},
i.e., that the user will not explicitly make space outside
span of the data to place the key. In this case, the user
should explicitly place the key as usual.
}
\seealso{
\code{\link{autoText}}, \code{\link{sparsestQuadrant}}, \code{\link{mergeList}}.}

\examples{
zkey <- list(x=-1,y=-1,
lines=list(type=c("o","l"), pch=c(18,0), lty=1, col=1:2),
border=1,
corner=c(1,1),
cex=1.5 )

old.mfrow <- par(mfrow=c(2,2))

## rectangular hole in cloud 
set.seed(100)
N <- 10024
y <- rnorm(N)
x <- seq(N)
xhole <- x > .4*N & x < .9*N
yhole <- y > -2 & y < 2
hole  <- xhole & yhole
y <- y[!hole]
x <- x[!hole]
plot(x, y, type="p")
autoKey(x, y, args=zkey)

## linear chirp (more room on left) 
y <- sin(.125*pi*100*((0:99)/100)^2)
plot(y, type="p")
autoKey(y, args=zkey)

## quadratic 
x <- seq(-5,5)
y <- x^2
plot(x,y,type="p")
autoKey(x, y, args=zkey)

## circle 
tt <- seq(0,2*pi,length=100)
x  <- sin(tt)
y  <- cos(tt)
plot(x, y, type="p")
autoKey(x, y, args=zkey)
par(old.mfrow)
}
\keyword{hplot}

