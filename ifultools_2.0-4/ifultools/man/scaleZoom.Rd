%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the scaleZoom function

\name{scaleZoom}
\alias{scaleZoom}
\title{Zoom in on a specified region of a data matrix}
\usage{scaleZoom(x, y, z=NULL, zoom=NULL, logxy="y",
    xy.linked=FALSE, xlab="x", ylab="y", log.base=2)}
\description{Returns the portion of a matrix based on the specified
ranges of abscissa (x) and ordinate (y) values.}

\arguments{
\item{x}{a numeric vector of abscissa values corresponding to the provided \code{z} matrix.}

\item{y}{a numeric vector of ordinate values corresponding to the provided \code{z} matrix.}
\item{log.base}{the base of the logarithm used to scale the output zoomed coordinates. See the help
for the \code{logxy} input for more details. Default: 2.}

\item{logxy}{a character string indicating containing either or both of the characters \code{"x"}
and \code{"y"}. If the character \code{"x"} is present in the string, the logarithm of the resulting
zoomed \code{x} coordinate is returned use base \code{log.base}. A similar logic follows for the
presence of character \code{"y"} in the string. For example, to take the log10 of both the zoomed x- and
y-coordinates, one would set \code{logxy="xy", log.base=10}. If logarithms are used, the \code{xlab} and
\code{ylab} strings are altered accordingly. Default: \code{"y"}.}

\item{xlab}{character string defining x-label. Default: \code{"x"}.}

\item{xy.linked}{a logical value. If \code{TRUE}, the algorithm used to
define the output (zoomed) x- and y-coordinates as follows:
\describe{
\item{1}{Let ix be the indices of \code{x} s.t. \code{x[ix]} is on the interval [\code{zoom[1]},\code{zoom[2]}].}

\item{2}{Redefine \code{x} and \code{y} as \code{x=x[ix]} and \code{y=y[ix]}.}

\item{3}{Let iy be the indices of \code{y} s.t. \code{y[iy]} is on the interval [\code{zoom[3]},\code{zoom[4]}].}

\item{4}{Redefine \code{x} and \code{y} as \code{x=x[iy]} and \code{y=y[iy]}.}}

Thus, both x- and y-zoom conditions must be met for both coordinates. Default: \code{FALSE}.}

\item{ylab}{character string defining x-label. Default: \code{"y"}.}

\item{z}{a matrix of numeric values whose x- and y-axis values are defined by the \code{x}
and \code{y} vectors, respectively. This argument may be \code{NULL} signifying that no \code{z} matrix
is present and (in this case) only the zoomed \code{x} and \code{y} coordinates will be returned
in the output. Default: \code{NULL}.}

\item{zoom}{a four-element numeric vector containing the ranges to zoom into the
data in \code{c(xmin, xmax, ymin, ymax)} format. Default: NULL (no zoom).}
}

\value{
a named list containing the following components:

\item{x}{The zoomed portion of the input \code{x} vector.}

\item{y}{The zoomed portion of the input \code{y} vector.}

\item{z}{The zoomed portion of the input \code{z} matrix, if \code{z} is not \code{NULL}.}

\item{ix}{The indices of the original \code{x} vector which span the range specified by \code{zoom[1:2]}.}

\item{iy}{The indices of the original \code{y} vector which span the range specified by \code{zoom[3:4]}.}

\item{xlab}{A character string defining an appropriate x-label for the zoomed x-coordinate output.}

\item{ylab}{A character string defining an appropriate y-label for the zoomed y-coordinate output.}
}
\seealso{
\code{\link[graphics]{image}}, \code{\link{stackPlot}}.}

\examples{
## create an image 
x <- y <- seq(-4*pi, 4*pi, len=50)
r <- sqrt(outer(x^2, y^2, "+"))
z <- cos(r^2)*exp(-r/6)
image(x,y,z,col=gray((0:32)/32))

## zoom in on one portion of that image and 
## re-display 
zoom <- scaleZoom(x,y,z, zoom=c(-5,5,-1,10), logxy="")
image(zoom$x, zoom$y, zoom$z, col=gray((0:32)/32))
}
\keyword{hplot}

