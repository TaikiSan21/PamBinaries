%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the linearFit function

\name{linearFit}
\alias{linearFit}
\title{Linear regression of a segmented time series}
\concept{segmentation, linear}
\usage{linearFit(x, y, fit=lmsreg, method="widest",
    n.fit=5, angle.tolerance=5, aspect=TRUE)}
\description{Segments the time series in approximately linear regions
as defined by the \code{linearSegmentation} function
and subsequently regressed the specified segment using
a selected regression scheme.}

\arguments{
\item{x}{the independent variable of the time series.}

\item{y}{the dependent variable of the time series.}
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

\item{fit}{a function representing the linear regression scheme to use in fitting
the resulting statistics (on a log-log scale). Supported functions are: \code{lm},
\code{lmsreg}, and \code{ltsreg}. See the on-line help documentation for each of these for
more information. Default: \code{lmsreg}.}

\item{method}{a character string used to define the criterion
for selecting one of the segments returned by the \code{linearSegmentation} function.
Choices are
\describe{
\item{\code{"first"}}{The first (leftmost) segment.}

\item{\code{"last"}}{The last (rightmost) segment.}

\item{\code{"widest"}}{The segment containing the most number of points.}

\item{\code{"all"}}{A union of all segments.}}

Default: \code{"widest"}.}

\item{n.fit}{an integer denoting the window size, not to exceed the
number of samples in the time series. Default: 5.}
}

\value{
the regression model.
}
\seealso{
\code{\link{logScale}}.}

\examples{
## obtain some data with approximately 
## piecewise-linear trends 
x <- seq(0,2*pi,length=100)
y <- sin(x)

## perform linear segmentation with aspect ratio 
## dilation using a 5 degree tolerance and 5 
## point windows. regress the widest of these 
## segments with the lm() function. 
linearFit(x, y, n.fit=5, angle.tolerance=5, aspect=TRUE,
    method="widest", fit=lm)
}
\keyword{utilities}

