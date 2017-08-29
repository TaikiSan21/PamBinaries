%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the splitplot function

\name{splitplot}
\alias{splitplot}
\title{Partitions plot space into rectangular grid}
\concept{utilities}
\usage{splitplot(nrows, ncols, i=1, new=as.logical(i > 1 && i <= nrows*ncols), gap=0.15)}
\description{Uses the \code{plt} argument of the \code{par} function to divide the space
according to the input grid.}

\arguments{
\item{nrows}{an integer defining the number of desired rows in the plot grid.}

\item{ncols}{an integer defining the number of desired column in the plot grid.}
\item{gap}{a numeric scalar used a nudge factor for the gap between plots in both the
x- and y-directions. Default: \code{0.15}.}

\item{i}{selects the \code{i}th plot of the current \code{nrow} by \code{ncol}
plot grid for the next plot. The upper left plot region is denoted by \code{i=1} and increments moving from
left to right, then top to bottom. Default: \code{1}.}

\item{new}{a logical flag. If \code{TRUE}, a new plot grid is established and previously defined plot regions
are erases. Default: \code{as.logical(i > 1)}.}
}

\value{
the original \code{plt} value of the \code{par} function prior to the call.
The user can use this to reset \code{plt} to its original values.
}
\note{The \code{plt} option of the \code{par} function is altered.}

\seealso{
\code{\link{stackPlot}}.}

\examples{
## establish a 2x2 plot grid and select the first 
## for plotting 
old.plt <- splitplot(2,2,1)

## plot the data 
for (i in seq(4)){
    if (i > 1)
       splitplot(2,2,i)

    plot(rnorm(100)*i)
    mtext(paste("i=", i, sep=""), side=3, line=0.5)
}

## return the original plot state of plt in par 
par(old.plt)
}
\keyword{hplot}

