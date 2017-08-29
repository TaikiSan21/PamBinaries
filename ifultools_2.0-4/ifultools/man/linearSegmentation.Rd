%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the linearSegmentation function

\name{linearSegmentation}
\alias{linearSegmentation}
\title{Piecewise linear segmentation of a time series}
\concept{segmentation, linear}
\usage{linearSegmentation(x, y, n.fit=5, angle.tolerance=5,
    aspect=TRUE, plot=FALSE, add=FALSE, ...)}
\description{Locates the change-points of time series based on a piecewise linear
segmentation algorithm. Given a window size (\code{n.fit}) and an angle tolerance
(\code{angle.tolerance}), the segmentation algorithm starts by finding the
slope of the first \code{n.fit} points of the series via least squares
regression. The window is slid over one point to the right, the points within
the new window are regressed, and the new slope is compared to the old slope.
If the change in slope exceeds the specified \code{angle.tolerance}, a change-point
is recorded as the rightmost point of the previous iteration's window. The
routine then picks up again starting at the point just to the right of the
change-point. If the change in slope does not exceed the specified
\code{angle.tolerance}, then the old slope is updated (in a running average sense),
and the algorithm continues as usual. The window is slid along the
until its rightmost point reaches the edge of the time series.}

\arguments{
\item{x}{the independent variable of the time series.}

\item{y}{the dependent variable of the time series.}
\item{...}{additional arguments sent unaltered to the \code{abline} function used to draw
the vertical divisions of the linear segmentation. Only
used if \code{plot=TRUE}.}

\item{add}{a logical value. If \code{TRUE}, the plot is added using the current plotting
scheme. Default: \code{FALSE}.}

\item{angle.tolerance}{the maximum angle in degrees that the running
average of the slopes in the current set of points
can change relative to the slope of the data calculated
in the most current (rightmost) window before a
change-point is recorded. Default: 5.}

\item{aspect}{a logical value. If \code{TRUE}, the aspect ratio
of the data (defined by \code{range(y) / range(x)}) is used to adjust
the \code{angle.tolerance}. Specifically, the new angle tolerance becomes
\code{angle.tolerance / aspect.ratio}. Using the
aspect ratio to dilate \code{angle.tolerance} allows the user
to specify the degree of variablity they wish to impose
on the data vis-a-vis a standard plot of the data, i.e. what
you would see if you issued \code{plot(xdata, ydata)}.
The idea is that when looking at such plots, one might
decide (for example) that a 5 degree variablity on the plot
would be acceptable. But if that range of \code{y} is vastly
different from that of \code{x}, then the true change in angle
from one section to the other will be much different than 5 degrees.
Thus, \code{aspect} can be used to compensate for aspect ratios far away from unity.
Default: \code{TRUE}.}

\item{n.fit}{an integer denoting the window size, not to exceed the
number of samples in the time series. Default: 5.}

\item{plot}{a logical value. If \code{TRUE}, a plot of the segmentation is displayed.
Default: \code{FALSE}.}
}

\value{
a vector containing the result.
}
\seealso{
\code{\link[stats]{lm}}.}

\examples{
## obtain some data with approximately 
## piecewise-linear trends 
x <- seq(0,2*pi,length=100)
y <- sin(x)

## perform linear segmentation with aspect ratio 
## dilation using a 5 degree tolerance and 5 
## point windows 
z <- linearSegmentation(x, y, n.fit=5, angle.tolerance=5, aspect=TRUE)

## plot the data and the estimated change-points 
plot(x, y)
abline(v=x[z], lty=2)
}
\keyword{data}
\keyword{manip}
\keyword{ts}

