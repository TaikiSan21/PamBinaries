%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the gridOverlay function

\name{gridOverlay}
\alias{gridOverlay}
\title{Add a grid over the existing plot}
\usage{gridOverlay(lty = "dotted", col = "gray", density = 3, ...)}
\description{Overlays a grid on the current plot at the ordinate and abscissa tick markers.}

\arguments{

\item{...}{additional arguments passed directly to the \code{par} function.}

\item{col}{line color ala \code{par} function. Default: \code{"gray"}.}

\item{density}{line density ala \code{par} function. Default: 3.}

\item{lty}{line type ala \code{par} function. Default: \code{"dotted"}.}
}
\seealso{
\code{\link[graphics]{par}}.}

\examples{
## plot the chaotic beam data and overlay a grid 
old.plt <- splitplot(1,1,1)
plot(rnorm(100))
gridOverlay()
par(old.plt)
}
\keyword{hplot}

