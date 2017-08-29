%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the solveODE function

\name{solveODE}
\alias{solveODE}
\title{Numerical integration of ODEs}
\concept{utilities}
\usage{solveODE(FUN, initial=NULL, step=0.01, n.sample=1000, n.transient=100, ...)}
\description{Solve a system of ordinary differential equations via a fourth order Runge-Kutta numerical
integration scheme.}

\arguments{
\item{FUN}{a function defining the sytem of ODEs. This function should have as an
input X, where X is a vector whose length is equal to the order of the ODEs. It should
return a value for each order (state) of the system.}
\item{...}{additional arguments sent directly to \code{FUN}.}

\item{initial}{a vector of initial values, one element for each state of the system defined by the ODEs.
By default, this value is \code{NULL}, in which case \code{FUN} is called with a \code{X=NULL}. Therefore,
\code{FUN} should be able to handle a \code{NULL} value input if you do not specify an initial condition. Default: \code{NULL}.}

\item{n.sample}{the number of desired samples for each state beyond that specified for the
transient ala \code{n.transient}. Default: 1000.}

\item{n.transient}{the numeber of transient points. These points are excluded from the
output. Default: 100.}

\item{step}{a numerical integration time step. Default: 0.1.}
}

\value{
a \code{data.frame} containing the estimated system response variables.
}
\seealso{
\code{\link[graphics]{par}}.}

\examples{
## estimate response of the chaotic Lorenz system 
"lorode" <- function(x, sigma = 10, r = 28, b = 8/3){
    c(sigma * (x[2] - x[1]), x[1] * (r - x[3]) - x[2],  - b * x[3]
       + x[1] * x[2])
}

z <- solveODE(lorode,  initial=c(0.1,0.3,1), n.transient=1500,
    n.sample=2000)
nms <- c("X","Y","Z")

## plot the results 
stackPlot(x=seq(150, by=0.1, length=2000), y=z,
    ylab=nms, main="Lorenz System in Chaos", xlab="Time")
}
\keyword{math}

