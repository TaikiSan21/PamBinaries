%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the ACVStoPACS function

\name{ACVStoPACS}
\alias{ACVStoPACS}
\title{Computes partial autocorrelations from autocovariances}
\concept{conversion}
\usage{ACVStoPACS(acvs)}
\description{Given an autocovariance sequence (ACVS) for a stationary process,
computes the corresponding partial autocorrelation sequence (PACS).}

\arguments{
\item{acvs}{the ACVS for lags 0, 1, ..., p, where p must be a positive integer.}

}

\value{
a vector with the PACS for lags 1, ..., p.
}
\details{
%
The PACS
(sometimes called the reflection coefficient sequence)
is computed from the ACVS using the Levinson-Durbin recursions.
Note that the autocorrelation sequence can be used as
input rather than the ACVS since the PACS does not
in fact depend on the variance of the process
(i.e., the ACVS at lag 0).
}
\references{
S-Plus 5 Guide to Statistics,Section 24.2.

D. Percival and A. Walden, \emph{Spectral Analysis for Physical Applications},
Cambridge University Press, 1993, Section 9.4.
}

\seealso{
\code{\link[stats]{ar.yw}}.}

\examples{
ACVStoPACS(c(3,2,1))
ACVStoPACS(c(1,0.9,0.81,0.9^3,0.9^4))
}
\keyword{utilities}

