%%  WARNING: This file was automatically generated from the associated 
%%  sapa_sdf.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the sapa_sdf.mid file in the project MID directory. Once the 
%%  sapa_sdf.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the SDF, as.matrix.SDF, plot.SDF, print.SDF functions

\name{SDF}
\alias{SDF}
\alias{as.matrix.SDF}
\alias{plot.SDF}
\alias{print.SDF}
\title{Nonparametric (cross) spectral density function estimation}
\concept{spectral density function estimation}
\usage{SDF(x, method="direct", taper.=NULL, window=NULL,
    n.taper=5, overlap=0.5, blocksize=NULL,
    single.sided=TRUE, sampling.interval=NULL,
    center=TRUE, recenter=FALSE, npad=2*numRows(x))}
\description{Estimate the process (cross) spectral density function via nonparametric
models.}

\arguments{
\item{x}{a vector or matrix containing uniformly-sampled real-valued time series.
If a \code{matrix}, each column should contain a different time series.}
\item{blocksize}{an integer representing the number
of points (width) of each block in the WOSA estimator scheme.
Default: \code{floor(N/4)} where \code{N} is the number of
samples in each series.}

\item{center}{a logical value. If \code{TRUE}, the mean of each
time series is recentered prior to estimating
the SDF. Default: \code{TRUE}.}

\item{method}{a character string denoting the method to use in estimating the SDF.
Choices are \code{"direct"}, \code{"lag window"}, \code{"wosa"} (Welch's Overlapped Segment Averaging),
\code{"multitaper"}. See \bold{DETAILS} for more information. Default: \code{"direct"}.}

\item{n.taper}{an integer defining the number of tapers
to use in a multitaper scheme. This value is overwritten if
the \code{taper} input is of class \code{taper}. Default: \code{5}.}

\item{npad}{an integer representing the total length of each
time series to analyze after padding with zeros. This argument
allows the user to control the spectral resolution of the SDF
estimates: the normalized frequency interval is
\eqn{\Delta f = 1 / \hbox{npad}}{deltaf=1/npad}.
This argument must be set such that
\eqn{\hbox{npad} > 2}{npad > 2}.
Default: \code{2*numRows(x)}.}

\item{overlap}{a numeric value on \eqn{[0,1]} denoting the fraction
of window overlap for the WOSA estimator. Default: \code{0.5}.}

\item{recenter}{a logical value. If \code{TRUE}, the mean of each
time series is recentered after (posssibly) tapering the series prior to estimating
the SDF. Default: \code{FALSE}.}

\item{sampling.interval}{a numeric value representing the interval
between samples in the input time series \code{x}. Default: \code{NULL}, which
serves as a flag to obtain the sampling interval via the \code{deltat}
function. If \code{x} is a list, the default sampling interval is \code{deltat(x[[1]])}.
If \code{x} is an atomic vector (ala \code{isVectorAtomic}), then the default
samplign interval is established ala \code{deltat(x)}. Finally, if the
input series is a matrix, the sampling interval of the first series (assumed to
be in the first column) is obtained ala \code{deltat(x[,1])}.}

\item{single.sided}{a logical value. If \code{TRUE}, a single-sided
SDF estimate is returned corresponding to the normalized frequency
range of \eqn{[0,1/2]}. Otherwise, a double-sided SDF estimate
corresponding to the normalized frequency interval \eqn{[-1/2,1/2]}
is returned. Default: \code{TRUE}.}

\item{taper.}{an object of class \code{taper} or a character string denoting the primary taper.
If an object of class \code{taper}, the length of the taper is checked to ensure
compatitbility with the input \code{x}.
See \bold{DETAILS} for more information. The default values are a function
of the \code{method} as follows:
\describe{
\item{direct}{normalized rectangular taper}

\item{lag window}{normalized Parzen window with a cutoff at \eqn{N/2} where
\eqn{N} is the length of the time series.}

\item{wosa}{normalized Hanning taper}

\item{multitaper}{normalized Hanning taper}}}

\item{window}{an object of class \code{taper} or a character string denoting the (secondary) window
for the lag window estimator.
If an object of class \code{taper}, the length of the taper is checked to ensure
compatitbility with the input \code{x}. See \bold{DETAILS} for more information.
Default: Normalized Hanning window.}
}

\value{
an object of class \code{SDF}.
}

\section{S3 METHODS}{
\describe{

\item{as.matrix}{converts the (cross-)SDF estimate(s) as a matrix. Optional arguments
are passed directly to the \code{matrix} function during the conversion.}

\item{plot}{plots the (cross-)SDF estimate(s). Optional arguments are:
\describe{
\item{xscale}{a character string defining the scaling to perform on the (common) frequency vector
of the SDF estimates. See the \code{scaleData} function for supported choices. Default: \code{"linear"}.}

\item{yscale}{a character string defining the scaling to perform on the SDF estimates.
See the \code{scaleData} function for supported choices. Default: \code{"linear"}.}

\item{type}{a single character defining the plot type (ala the \code{par} function) of the
SDF plots. Default: \code{ifelse(numRows(x) > 100, "l", "h")}.}

\item{xlab}{a character string representing the x-axis label. Default: \code{"FREQUENCY (Hz)"}.}

\item{ylab}{a (vector of) character string(s), one per (cross-)SDF estimate,
representing the y-axis label(s). Default: in the multivariate case, the strings
\code{"Sij"} are used for the y-axis labels, where i and j are the indices of the
different variables. For example, if the user supplies a 2-column matrix for \code{x},
the labels \code{"S11"}, \code{"S12"}, and \code{"S22"} are used to label the y-axes of the corresponding
(cross-)SDF plots. In the univariate case, the default string \code{"SDF"} prepended with a string
describing the type of SDF performed (such as \code{"Multitaper"}) is used to label the y-axis.}

\item{plot.mean}{a logical value. If \code{TRUE}, the SDF value at normalized frequency \eqn{f=0}
is plotted for each SDF. This frequency is associated with the sample mean
of the corresponding time series. A relatively large mean value dominates
the spectral patterns in a plot and thus the corresponding frequency is typically not plotted.
Default: \code{!attr(x,"center")}.}

\item{n.plot}{an integer defining the maximum number of SDF plots to place onto a single graph.
Default: \code{3}.}

\item{FUN}{a post processing function to apply to the SDF values prior to plotting. Supported
functions are \code{Mod}, \code{Im}, \code{Re} and \code{Arg}. See each of these functions for details.
If the SDF is purely real (no cross-SDF is calculated), this argument is coerced to the \code{Mod} function.
Default: \code{Mod}.}

\item{add}{A logical value. If \code{TRUE}, the plot is added using the current
\code{par()} layout. Otherwise a new plot is produced. Default: \code{FALSE}.}

\item{...}{additional plot parameters passed directly to the \code{genPlot} function used
to plot the SDF estimates.}}}

\item{print}{prints the object. Available options are:
\describe{
\item{justify}{text justification ala \code{prettPrintList}. Default: \code{"left"}.}

\item{sep}{header separator ala \code{prettyPrintList}. Default: \code{":"}.}

\item{...}{Additional print arguments sent directly to the \code{prettyPrintList} function.}}}

}
}

\details{
%
Let \eqn{X_t}{x(t)} be a uniformly sampled real-valued time series of length \eqn{N},
Let an estimate of the process spectral density function be denoted as
\eqn{\hat{S}_X(f)}{S(f)} where \eqn{f} are
frequencies on the interval
\eqn{[-1/(2\Delta t),1/(2\Delta t)]}{-1/(2*deltat),1/(2*deltat)} where
\eqn{\Delta t}{deltat} is the sampling interval. The supported SDF
estimators are:

\describe{
\item{direct}{The direct SDF estimator is defined as
\eqn{\hat{S}_X^{(d)}(f) = | \sum_{t=0}^{N-1} h_t X_t e^{-i2\pi f t}|^2}{S(f)=|sum[t=0,...,N-1]{h(t)*x(t)*exp(-i*2*pi*f*t)}|^2},
where \eqn{\{h_t\}}{h(t)} is a data taper normalized such that
\eqn{\sum_{t=0}^{N-1} h_t^2 = 1}{sum[t=0,...,N-1]{h(t)^2} = 1}. If
\eqn{h_t=1/\sqrt{N}}{h(t)=1/sqrt(N)} then we obtain the definition
of the periodogram
\eqn{\hat{S}_X^{(p)}(f) = \frac{1}{N} | \sum_{t=0}^{N-1} X_t e^{-i2\pi f t}|^2}{S(f)=(1/N) * |sum[t=0,...,N-1]{x(t)*exp(-i*2*pi*f*t)}|^2}.
See the \code{taper} function for more details on supported window types.}

\item{lag window}{The lag window SDF estimator is defined as
\eqn{\hat{S}_X^{(lw)}(f) = \sum_{\tau=-(N-1)}^{N-1} w_\tau
\hat{s}_{X,\tau}^{(d)} e^{-i2\pi f \tau}}{S(f)=sum[k=-(N-1),...,(N-1)]{w(k)*s(k)*exp(-i*2*pi*f*k)}|^2},
where \eqn{\hat{s}_{X,\tau}^{(d)}}{s(k)} is the autocovariance
sequence estimator corresponding to some
direct spectral estimator (often the periodogram)
and \eqn{w_\tau}{w(k)} is a lag window (popular choices
are the Parzen, Papoulis, and Daniell windows). See the
\code{taper} function for more details.}

\item{wosa}{Welch's Overlapped Segment Averaging SDF estimator is
defined as

\deqn{
\hat S^{(wosa)} =
{1\over N_B} \sum_{j=0}^{N_B-1}
\hat S^{(d)}_{jN_O} (f)
}{S(f)=(1/Nb)*sum[j=0,...,Nb-1]{S(j*No,f)}}
where
\deqn{
\hat S^{(d)}_{l}(f)
\equiv
\left|
\sum_{t=0}^{N_S-1}
h_t X_{t+l} e^{-i2\pi ft}
\right|^2,
\enskip
0 \le l \le N - N_S;
}{S(l,f) =|sum[t=0,...,Ns-1]{h(t)*x(t+l)*exp(-i*2*pi*f*t)}|^2}
Here, \eqn{N_O}{No} is a positive integer
that controls how much overlap there is between segments
and that must satisfy both
\eqn{N_O \le N_S}{No <= Ns} and \eqn{N_O(N_B-1) = N-N_S}{No * (Nb - 1) = N - Ns},
while \eqn{\{ h_t \}}{h(t)} is a data taper appropriate
for a series of length \eqn{N_S}{Ns}
(i.e., \eqn{\sum_{t=0}^{N_S-1} h_t^2 = 1}{sum[t=0,...,Ns-1]{h_t^2} = 1}).}

\item{multitaper}{A multitaper spectral estimator is given by

\deqn{\hat S^{(mt)}_X(f)=
{1\over K}
\sum_{k=0}^{K-1}
\left| \sum_{t=0}^{N-1} h_{k,t} X_t e^{-i2\pi ft} \right|^2,
}{S(f) = (1/K) * sum[k=0,...,K-1] S(k,f)}
where
\eqn{S(k,f) = {|\sum_{t=0}^{N-1} h_{k,t} X_t  \exp(-i 2 \pi f t)|}^2}{S(k,f) = |sum[t=0,...,N-1]{h(k,t) * X(t) * exp(-i*2*pi*f*t)}|^2}
and \eqn{\{ h_{k,t} \}$, $k=0,\ldots,K-1}{h(k,t) for k=0,...,K-1},
is a set of \eqn{K} orthonormal data tapers.
\deqn{\sum_{t=0}^{N-1} h_{k,t} h_{k',t}
=
\left\{
\begin{array}{ll}
1,& \mbox{if $k=k'$;}\\
0,& \mbox{otherwise}
\end{array}
\right.
}{See reference(s) for further details.}
Popular choices for multitapers include sinusoidal tapers and
discrete prolate spheroidal sequences (DPSS). See the
\code{taper} function for more details.}}


\bold{Cross spectral density function estimation:}
If the input \code{x} is a matrix, where each column contains
a different time series, then the results are returned in
a matrix whose columns correspond to all possible unique combinations
of cross-SDF estimates. For example, if \code{x} has three columns,
then the output will be a matrix whose columns are
\eqn{\{S_{11},S_{12},S_{13},S_{22},S_{23},S_{33}\}}{{S11, S12, S13, S22, S23, S33}}
where \eqn{S_{ij}}{Sij} is the cross-SDF estimate of the \code{i}th
and \code{j}th column of \code{x}. All cross-spectral density
function estimates are returned as complex-valued series to maintain
the phase relationships between components. For all \eqn{S_{ij}}{Sij}
where \eqn{i=j}, however, the imaginary portions will be zero (up to a
numerical noise limit).
}
\references{
Percival, Donald B. and Constantine, William L. B. (2005)
``Exact Simulation of Gaussian Time Series
from Nonparametric Spectral Estimates
with Application to Bootstrapping",
\emph{Journal of Computational and Graphical Statistics},
accepted for publication.

D.B. Percival and A. Walden (1993),
\emph{Spectral Analysis for Physical Applications: Multitaper
and Conventional Univariate Techniques},
Cambridge University Press, Cambridge, UK.
}

\seealso{
\code{\link{taper}}, \code{\link{ACVS}}.}

\examples{
## calculate various SDF estimates for the 
## sunspots series. remove mean component for a 
## better comparison. 

require(ifultools)
data <- as.numeric(sunspots)
methods <- c("direct","wosa","multitaper",
    "lag window")

S <- lapply(methods, function(x, data) SDF(data, method=x), data)
x <- attr(S[[1]], "frequency")[-1]
y <- lapply(S,function(x) decibel(as.vector(x)[-1]))
names(y) <- methods

## create a stack plot of the data 
stackPlot(x, y, col=1:4)

## calculate the cross-spectrum of the same 
## series: all spectra should be the same in 
## this case 
SDF(cbind(data,data), method="lag")

## calculate the SDF using npad=31 
SDF(data, npad=31, method="multitaper")
}
\keyword{univar}

