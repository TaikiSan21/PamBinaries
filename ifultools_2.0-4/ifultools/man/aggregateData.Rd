%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the aggregateData function

\name{aggregateData}
\alias{aggregateData}
\title{Time series aggregation}
\concept{data aggregation}
\usage{aggregateData(x, by, FUN, moving=FALSE, ...)}
\description{Partitions data into blocks and and applies the specifed function to each of the bocks.}

\arguments{
\item{x}{a numeric vector.}

\item{by}{.}

\item{FUN}{a scalar function to compute the summary statistics which can be applied to all data subsets.}
\item{...}{additional arguments to pass to FUN.}

\item{moving}{either FALSE to do standard aggregation, or a positive integer N to perform a moving aggregation (normally used for a moving average) over N samples.}
}

\value{
the aggregated series.
}
\seealso{
\code{\link[stats]{aggregate}}.}

\examples{
## Group a simple series into blocks containing 8 
## elements, and take the mean of each block. 
## Each block is lagged by 3 elements 
aggregateData(1:30, by=3, FUN=mean, moving=8)
}
\keyword{utilities}

