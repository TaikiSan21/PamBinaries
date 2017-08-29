%%  WARNING: This file was automatically generated from the associated 
%%  it_mutils.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_mutils.mid file in the project MID directory. Once the 
%%  it_mutils.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mutilsSDF function

\name{mutilsSDF}
\alias{mutilsSDF}
\title{SDF evaluation over a discrete uniform grid of frequencies}
\concept{utilities}
\usage{mutilsSDF(sdf=NULL, sdfargs=NULL, n.freq=1024, sampling.interval=1)}
\description{MUTILS is a C library used for certain \code{itCall} functions.
This function evaluates an input SDF function over frequencies on
the interval [0, Nyquist], i.e., normalized frequencies [0, 1/2].}

\arguments{

\item{n.freq}{a positive integer (greater than one) defining the number of frequencies to use in evaluating the SDF. The
frequencies are uniformly distributed over the interval [0, Nyquist] ala f=[0, 1/P , 2/P, 3/P, ..., (n.freq-1)/P]
where P=2*(n.freq-1)/sampling.interval. Default: \code{1024}.}

\item{sampling.interval}{a positive numeric scalar representing the sampling interval of the time series
associated with the input SDF. This argument is used only develop the set of frequencies over which the input
SDF is evaluated (see documentation for \code{n.freq} argument for more details). Default: \code{1}.}

\item{sdf}{the input SDF. If \code{sdf=NULL}, a negative numeric scalar is returned and serves as a flag the MUTILS C code to ignore
relevant SDF calculations. Otherwise, this input must be a function whose first argument is \code{f}
(representing frequency). At a minimum, the SDF must be defined over frequencies [0, Nyquist]
where \code{Nyquist=1/(2*sampling.interval)}. Additional arguments that are needed to calculate the SDF should
be passed via the \code{sdfargs} parameter. Default: \code{NULL}.}

\item{sdfargs}{a list of named additional SDF arguments. For example, if the input SDF
requires additional arguments (other than \code{f}): \code{y} and \code{z}, then specify this input
ala \code{sdfargs=list(y=1, z=pi)}. Default: \code{NULL} (no additional inputs).}
}

\value{
a numeric vector containing the values of the input SDF evaluated over frequencies uniformly distributed on the interval [0, Nyquist].
The corresponding frequencies are assigned to the output object as the named attribute \code{"frequency"}.
If \code{sdf=NULL}, a negative numeric scalar is returned and serves as a flag to MUTILS C code that the SDF is missing or undefined.
}
\seealso{
\code{\link{mutilsDistanceMetric}}, \code{\link{mutilsFilterType}}, \code{\link{mutilsFilterTypeContinuous}}, \code{\link{mutilsTransformPeakType}}, \code{\link{mutilsTransformType}}, \code{\link{mutilsWSShrinkageFunction}}, \code{\link{mutilsWSThresholdFunction}}.}

\examples{
## create a faux SDF 
S <- function(f, phi)  1/(1 + phi^2 - 2*phi*cos(2*pi*f))

## specify additional input arguments needed to 
## evaluate the SDF 
sdfargs <- list(phi=0.9)

## evaluate the SDF over frequencies on the 
## interval [0, 1/2] 
Sx <- mutilsSDF(sdf=S, sdfargs=sdfargs)

## plot the result 
f <- attr(Sx,"frequency")
plot(f, Sx, type="l")
}
\keyword{IO}

