%%  WARNING: This file was automatically generated from the associated 
%%  it_plot.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_plot.mid file in the project MID directory. Once the 
%%  it_plot.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the stackPlot function

\name{stackPlot}
\alias{stackPlot}
\title{Stack plot}
\concept{utilities}
\usage{stackPlot(x, y = NULL, xlty = NULL, bty = "n", lty = 1, col = 1:8,
    lwd = 1, rescale = TRUE, add = FALSE, cex = 1, xaxs = "r", xpd = TRUE, 
    yaxis = list(add = TRUE, ndigit = 3,
    col = 1:8, lty = 1, lwd = 3, side = "left", cex = 1),
    xlab = list(text = "", cex = 1, srt = 0, col = 1), ylab = list(text = NULL,
    cex = 1, srt = 0, col = 1:8, side = "right"), main = list(text = "",
    cex = 1, srt = 0, col = 1, adj = 0.5), ylim = NULL)}
\description{Plots input data as a stack of plots.}

\arguments{
\item{x}{a vector of numeric values corresponding to a common abscissa (x-axis) for
all ordinate (y-axis) values.}

\item{y}{a numeric vector, matrix, data.frame, or rectangular list
containing the ordinate (y-axis) values.}
\item{add}{a logical value. If \code{TRUE}, the plot is added using the current
\code{par} layout. Otherwise a new plot is produced. Default: \code{FALSE}.}

\item{bty}{a character string defining the box type (ala \code{par}'s \code{bty} parameter)
for each plot. Default: \code{"n"} (no boxes).}

\item{cex}{a numeric value defining the character expansion for the plot labels (ala \code{par}'s \code{cex} parameter).
Default: 1.}

\item{col}{an integer or vector of integers denoting the color of each
plotted series (ala \code{par}'s \code{col} parameter). This
input can be a vector, one color for each series plotted. If the length
of this vector is less than the number of series, then \code{col} is repeated as many
times as necessary. Default: \code{1:8}.}

\item{lty}{the line type for each plot (ala \code{par}'s \code{lty} parameter). This
input can be a vector, one line type for each series plotted. If the length
of this vector is less than the number of series, then \code{lty=1} is used for those
plots where lty is not directly specified. Default: 1.}

\item{lwd}{an integer or vector of integers denoting the line width of each
plotted series (ala \code{par}'s \code{lwd} parameter). This
input can be a vector, one line width for each series plotted. If the length
of this vector is less than the number of series, then \code{lwd} is repeated as many
times as necessary. Default: 1.}

\item{main}{the main label as defined by a list of the following named objects:
\describe{
\item{text}{A character string defining the label. Default: \code{""} (no label).}

\item{cex}{A character expansion value used to scale the label. Default: 1.}

\item{srt}{A numerical value specifying the rotation of the label in degrees. Default: 0.}

\item{col}{An integer defining the color of the label. Default: 1.}

\item{adj}{A numeric value on [0,1] defining the justification of the label relative to
the entire width of the plot window. The value 0, 0.5, and 1 represent left, center,
and right text alignment, respectively. Default: 0.5.}}

If a partial list of the above named objects is supplied, those objects
are merged with the default list defined above.}

\item{rescale}{a logical value. If \code{TRUE}, the data in each plot is scaled so that
visually the height of each subplot is approximately the same (the y-axis labels denote
the unscaled/original values. Default: \code{TRUE}.}

\item{xaxs}{a character string defining the style of the x-axis
interval calculation (ala \code{par}'s \code{xaxs} parameter). Default: \code{"e"} (extended axes).}

\item{xlab}{the x-axis label as defined by a list of the following named objects:
\describe{
\item{text}{A character string defining the label. Default: \code{""} (no label).}

\item{cex}{A character expansion value used to scale the label. Default: 1.}

\item{srt}{A numerical value specifying the rotation of the label in degrees. Default: 0.}

\item{col}{An integer defining the color of the label. Default: 1.}}

If a partial list of the above named objects is supplied, those objects
are merged with the default list defined above.}

\item{xlty}{the line type (ala \code{par}'s \code{lty} parameter) of
horizontal lines used to separate stacked plots. Default: \code{NULL} (no separator lines).}

\item{xpd}{a character string defining the style of the plot clipping
(ala \code{par}'s \code{xpd} parameter). Default: \code{TRUE}.}

\item{yaxis}{the y-axis style as defined by a list of the following named objects:
\describe{
\item{add}{A logical value. If \code{TRUE}, a y-axis is drawn for each plot. Default: \code{TRUE}.}

\item{ndigit}{An integer defining the number of digits to use in labeling y-axis data ranges. Default: 3.}

\item{lwd}{An integer defining the line width of the y-axis. Default: 3.}

\item{side}{A character string (either \code{"left"} or \code{"right"} denoting
the side to place the y-axis. Default: \code{"left"}.}

\item{cex}{A character expansion value used to scale the labels on the y-axis. Default: 1.}}

If a partial list of the above named objects is supplied, those objects
are merged with the default list defined above.}

\item{ylab}{the y-axis label(s) as defined by a list of the following named objects:
\describe{
\item{text}{A (vector of) character string(s) defining the label(s).
Default: the names of the y-axis labels as labeled in the primary \code{y} input.}

\item{cex}{A character expansion value used to scale the label(s). Default: 1.}

\item{srt}{A numerical value specifying the rotation of the label(s) in degrees. Default: 0.}

\item{col}{A (vector of) integer(s) defining the color of the label(s). Default: \code{1:8}.}

\item{side}{A character string (either \code{"left"} or \code{"right"} denoting
the side to place the y-axis label(s). Default: \code{"right"}.}}

If a partial list of the above named objects is supplied, those objects
are merged with the default list defined above.}

\item{ylim}{a two-element numeric vector containing the y-axis range of each series
to plot. The range of the specified \code{ylim} values must span that of all the series,
otherwise an error is returned. Default: \code{NULL} (y-axis limits are calculated automatically).}
}

\value{
no output.
}
\note{A stack plot is produced.}

\seealso{
\code{\link{splitplot}}.}

\examples{
## stack-plot the sunspots series and a random 
## walk series 
set.seed(100)
ix <- seq(2048)
stackPlot(x=ix,
    y=data.frame(sunspots[ix], cumsum(rnorm(length(ix)))),
    xlty=2, ylab=list(text=c("sunspots","walk")))
}
\keyword{hplot}

