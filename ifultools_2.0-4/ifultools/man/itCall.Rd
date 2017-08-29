\name{itCall}
\alias{itCall}
\title{
Thin itCall Wrapper to IFULTOOLS Symbols
}
\description{
Thin itCall Wrapper to IFULTOOLS Symbols
}
\usage{
itCall(symbol, ...)
}
\arguments{
  \item{symbol}{character scalar defining symbol to call in DLL}
  \item{\dots}{arguments to underlying C code}
}
\details{
Foreign function calls are no longer allowed in CRAN. This function serves as a thin wrapper to avoid 
associated R CMD check issues when building packages that depend on IFULTOOLS.
}

\value{output of the \code{itCall}}

\seealso{
\code{\link{itCall}}.
}

\examples{
\dontrun{
itCall("RS_fractal_filter_nonlinear_local_projection",
    as.matrix(x),
    as.integer(dimension),
    as.integer(tlag),
    as.integer(n.neighbor),
    max.distance,
    mutilsDistanceMetric(metric),
    as.integer(noise.dimension),
    corr.curve)
}
}

\keyword{utilities}
