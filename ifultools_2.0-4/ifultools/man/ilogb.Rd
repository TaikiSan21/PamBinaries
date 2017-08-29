%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the ilogb function

\name{ilogb}
\alias{ilogb}
\title{Integer truncation of logb output}
\concept{utilities}
\usage{ilogb(x, base=2, eps=.Machine$double.eps * 1e9)}
\description{Calculates \code{logb(x, base=base) + eps} where \code{eps} is a user-supplied (small) value and
converts the result ala \code{as.integer}. The purpose of this function is to avoid
\code{logb} variations between platforms. As an example, it is conceivable that \code{logb(32, base=2)}
on one platform could return a value of \code{5} while on another platform it could return a slightly smaller value,
e.g., \code{4.999999999999999}. The difference is due to different compiler and platform-specific implementations of
the S-PLUS or R languages. If \code{logb} output is truncated (as is the case with \code{as.integer(logb(x))}) then the
result may be unexpected, leading to intractable errors. By adding an appropriate and small value to \code{logb(x)},
subsequent truncation results in consistent value(s) across platforms.}

\arguments{
\item{x}{a numeric or complex vector.}
\item{base}{a positive or complex number: the base with respect to which logarithms are computed. Default: \code{2}.}

\item{eps}{a small numeric value that is added to the \code{logb} result prior to truncation via \code{as.integer}.
Default: \code{.Machine\$double.eps * 1e9}.}
}

\value{
a numeric scalar or vector containing the result.
}
\seealso{
\code{\link{logScale}}.}

\examples{
## should return 3 
ilogb(8 - .Machine$double.eps, base=2)
}
\keyword{utilities}

