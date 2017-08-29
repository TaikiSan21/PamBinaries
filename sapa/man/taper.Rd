%%  WARNING: This file was automatically generated from the associated 
%%  sapa_taper.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the sapa_taper.mid file in the project MID directory. Once the 
%%  sapa_taper.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the taper, as.matrix.taper, plot.taper, print.taper functions

\name{taper}
\alias{taper}
\alias{as.matrix.taper}
\alias{plot.taper}
\alias{print.taper}
\title{Oracle function for obtaining a particular taper/window}
\concept{taperswindows}
\usage{taper(type="rectangle", n.sample=100, n.taper=NULL,
    sigma=0.3, beta=4*pi*(n.sample-1)/n.sample, cutoff=floor(n.sample/2),
    sidelobedB=80, roughness=n.sample/2, flatness=0.3,
    bandwidth=4, normalize=TRUE)}
\description{Develop signal processing tapers or windows.}

\arguments{

\item{bandwidth}{bandwidth for DPSS tapers.
See \code{Details} for more information. Default: 4.}

\item{beta}{kaiser window shape factor (must be positive or zero).
See \code{Details} for more information. Default: \code{4*pi*(n.sample-1)/n.sample}.}

\item{cutoff}{parzen or Papoulis window cutoff (must be greater than unity).
See \code{Details} for more information. Default: \code{floor(n.sample/2)}.}

\item{flatness}{raised cosine taper flatness fraction (must be on [0,1]).
See \code{Details} for more information. Default: 0.3.}

\item{n.sample}{an integer denoting the number of samples.
Default: 1000.}

\item{n.taper}{an integer defining the multitaper order (number of orthogonal tapers) to use
in a multitaper scheme. The taper order directly impacts the quality of
the SDF estimate. Low taper orders are usually associated with SDF estimates with low bias and high variance,
while high taper orders attenuate the variance of the estimate at the risk of incurring a large bias.
This tradeoff between bias and variance is unavoidable but taper order allows you to tune the SDF
to meet the needs of your application. Studies show that a multitaper order of 5 typically provides a good balance
with reasonably low bias and variance properties (see the references for more details). Default: \code{NULL},
which serves as a flag to set the default taper order depending on the type of taper chosen for the analysis.
If \code{sine} or \code{dpss} multitapers are chosen, the default taper order is 5, otherwise is set to unity.}

\item{normalize}{a logical value. If \code{TRUE},
the taper is normalized to have unit energy. Default: \code{TRUE}.}

\item{roughness}{daniell windown roughness factor (must be positive).
See \code{Details} for more information. Default: \code{n.sample/2}.}

\item{sidelobedB}{chebyshev sidelobedB bandwidth in decibels
(must be positive).
See \code{Details} for more information. Default: 80.}

\item{sigma}{standard deviation for Guassian taper. Default: 0.3.}

\item{type}{a character string denoting the type of taper to create.
Supported types are \code{"rectangle"}, \code{"triangle"}, \code{"raised cosine"},
\code{"hanning"}, \code{"hamming"}, \code{"blackman"}, \code{"nuttall"}, \code{"gaussian"},
\code{"kaiser"}, \code{"chebyshev"}, \code{"born jordan"}, \code{"sine"},
\code{"parzen"}, \code{"papoulis"}, \code{"daniell"}, and \code{"dpss"}.
See \code{Details} for more information. Default: \code{"rectangle"}.}
}

\value{
an object of class \code{taper}.
}

\section{S3 METHODS}{
\describe{

\item{as.matrix}{converts output to a matrix.}

\item{plot}{plots the output. Optional arguments are:
\describe{
\item{ylab}{Character string denoting the y-axis label for
the plot. Default: \code{upperCase(attr(x,"type"))}.}

\item{type}{Line type (same as the \code{type} argument of the \code{par} function). Default: "l".}

\item{...}{Additional plot arguments (set internally by the \code{par} function).}}}

\item{print}{prints a summary of the output object.}

}
}

\details{
%
Let \eqn{w(\cdot)}{w()} and \eqn{h(t)} for t = 0, ..., N-1
be a \emph{lag window} and \emph{taper}, respectively.
The following lag window or taper types are supported.

\describe{
\item{rectangular}{A rectangular taper is defined as
\eqn{h_t=1}{h(t) = 1}.}

\item{triangle}{A triangular taper is defined as
\eqn{h_t=h_{N-t-1}=2(t+1)/(N+1)}{h(t)=h(N-t-1)=2*(t+1)/(N+1)}
for \eqn{t < M}{t < M} where \eqn{M = \lfloor N / 2 \rfloor}{M=floor(N/2)}
and \eqn{h_M=1}{h(M) = 1} if \eqn{N} is evenly divisible by 2.}

\item{raised cosine}{A raised cosine is a symmetric
taper with a flat mid-plateau. Let \eqn{p \in [0,1]}{p on [0,1]}
be the fraction of the length of the taper that is flat,
\eqn{M=\lfloor pN \rfloor}{M=floor(p*N)}, and
\eqn{\beta=2\pi/(M+1)}{B=2*pi/(M+1)}.
A raised cosine taper is defined as
\deqn{
h_t=h_{N-t-1}= 0.5 (1 - \cos(\beta(t+1))) \enskip\hbox{for}\enskip
0 \le t < \lfloor M/2 \rfloor
}{h(t) = h(N-t-1) = 0.5 * (1 - cos(B*(t+1))) for 0 <= t < floor(M/2) and}
\deqn{
h_t = 1 \enskip\hbox{for}\enskip \lfloor M/2 \rfloor \le t < N - \lfloor M/2 \rfloor.
}{h(t) = 1 for floor(M/2) <= t < N - floor(M/2).}}

\item{hanning}{Let \eqn{\beta=2\pi/(N+1)}{B=2*pi/(N+1)}.
A Hanning taper is defined as
\eqn{h_t=0.5 (1 - \cos(\beta(t+1)))}{h(t) = 0.5 * (1 - cos(B*(t+1)))}.}

\item{hamming}{Let \eqn{\beta=2\pi/(N-1)}{B=2*pi/(N-1)}.
A Hamming taper is defined as
\eqn{h_t=0.54 - 0.46\cos(\beta t)}{h(t) = 0.54 - 0.46 * cos(B*t)}.}

\item{blackman}{Let \eqn{\beta=2\pi/(N+1)}{B=2*pi/(N+1)}.
A Blackman taper is defined as
\eqn{h_t=0.42 - 0.5 \cos( \beta ( t + 1 ) ) +
0.08 \cos( 2 \beta ( t + 1 ) )}{h(t) = 0.42 - 0.5 * cos( B * ( t + 1 ) ) +
0.08 * cos( 2 * B * ( t + 1 ) )}.}

\item{nuttall}{Let \eqn{\beta=2\pi/(N-1)}{B=2*pi/(N-1)}.
A Nuttall taper is defined as
\eqn{h_t=0.3635819 - 0.4891775 \cos( \beta t ) + 0.1365995 \cos( 2\beta t )
- 0.0106411 \cos( 3\beta t )}{h(t) = 0.3635819 - 0.4891775 * cos( B*t ) + 0.1365995 * cos( 2*B*t ) -
0.0106411 * cos( 3*B*t )}}

\item{gaussian}{Let \eqn{\sigma}{S} be the standard deviation
of a Gaussian distribution. Let \eqn{\beta(t)=2\sigma(0.5 - t/(N-1))}{B(t)=2*S*(0.5 - t/(N-1))}.
A Gaussian taper is defined as
\deqn{h_t = h_{N-t-1} = e^{-B^2/2} \enskip\hbox{for}\enskip 0 \le t < \lfloor N/2 \rfloor}{h(t) = h(N-t-1) = exp(-B*B/2) for 0 <= t < floor(N/2)}.
\deqn{h_{N/2}=1\enskip\hbox{if N is evenly divisible by 2}}{h(N/2)=1 if N is evenly divisible by 2}.}

\item{kaiser}{Let \eqn{V_t = (2t-1-N)/N}{V(t) = (2t-1-N)/N} and \eqn{I_0(\cdot)}{I0()} be the zeroth-order modified
Bessel function of the first kind.  Given the shape factor \eqn{\beta > 0}{beta > 0},
a Kaiser taper is defined as
\eqn{h_t = I_0(\beta \sqrt{1-V_t^2})/I_0(\beta)}{h(t) = I0(beta*sqrt(1-V(t)\eqn{\mbox{\textasciicircum}}{^}2))/I0(beta)}.}

\item{chebyshev}{The Dolph-Chebyshev taper is a function of
both the desired length \eqn{N} and the desired sidelobe level
(our routine accepts a sidelobe attenuation factor expressed
in decibels). See the Mitra reference for more details.}

\item{born jordan}{Let \eqn{M=(N-1)/2}.
A Born-Jordan taper is defined as
\eqn{h_t==h_{N-t-1} = 1 / ( M - t + 1 )}{h(t) = h(N-t-1) = 1 / ( M - t + 1 )}.}

\item{sine}{Sine multitapers are defined as
\deqn{
h_{k,t} =
\left({2\over{N+1}}\right)^{1/2}
\sin\,\left( {{(k+1)\pi (t+1)} \over {N+1}}\right),
}{h(k,t) = sqrt(2/(N+1)) * sin((k+1)*pi*(t+1)/(N+1))}
for t = 0, ..., N-1 and k = 0, ..., .. This simple equation
defines a good approximation to the discrete prolate
spheroidal sequences (DPSS) used in multitaper SDF
estimation schemes.}

\item{parzen}{A Parzen lag window is defined as
\deqn{
w_{\tau,m}
=
\left\{
\begin{array}{lll}
{1 - 6(t/m)^2 + 6(|t|/m)^3},& \mbox{$|t| \le m/2$;}\\
{2(1 - |t|/m)^3},& \mbox{$m/2 < |t| \le m/2$;}\\
0,& \mbox{otherwise}.
\end{array}
\right.
}{w(t,m)= {1 - 6(t/m)^2 + 6(|t|/m)^3} for |t| <= m/2;
{2(1 - |t|/m)^3} for m/2 < |t| <= m; and 0 otherwise}
for \eqn{-(N-1) \le \tau \le (N - 1)}{-(N-1) <= t <= (N - 1)}.
The variable \eqn{m} is referred to as the \emph{cutoff} since
all values beyond that point are zero.}

\item{papoulis}{A Papoulis lag window is defined as
\deqn{
w_{\tau,m}
=
\left\{
\begin{array}{ll}
{{1 \over \pi} |\sin(\pi\tau/m)| +(1-|\tau|/m)\cos(\pi\tau/m)},& \mbox{$|\tau| < m$;}\\
0,& \mbox{$|\tau| \ge m$}
\end{array}
\right.
}{w(t,m)= {(1/pi) * |sin(pi*t/m)| + (1 - |t|/m) * cos(pi*t/m)} for |t| < m;
and 0 otherwise}
for \eqn{-(N-1) \le \tau \le (N - 1)}{-(N-1) <= t <= (N - 1)}.
The variable \eqn{m} is referred to as the \emph{cutoff} since
all values beyond that point are zero.}

\item{daniell}{A Daniell lag window is defined as
\deqn{
w_{\tau,m}
=
\left\{
\begin{array}{ll}
{\sin(\pi\tau/m) \over \pi\tau/m},& \mbox{$|\tau| < N$;}\\
0,& \mbox{$|\tau| \ge N$}
\end{array}
\right.
}{w(t,m)= {sin(pi*t/m)/(pi*t/m)} for |t| < N;
and 0 for |t| >= N}
for \eqn{-(N-1) \le \tau \le (N - 1)}{-(N-1) <= t <= (N - 1)}.
The variable \eqn{m} is referred to as the \emph{roughness} factor,
since, in the context of spectral density function (SDF) estimation,
it controls the degree of averaging that is performed
on the preliminary direct SDF estimate. The smaller the roughness,
the greater the amount of smoothing.}

\item{dpss}{Discrete prolate spheroidal sequences are (typically) used
for multitaper spectral density function estimation. The first order
DPSS can be defined (to a good approximation) as
\deqn{h_{t,0} = C \times I_0\bigl(\widetilde{W} \sqrt{1 - (1 - g_t)^2}\bigr)
/I_0(\widetilde{W} )
}{h(t,k) = C * I0( W' * sqrt(1 - (1 - g(t))^2) / I0(W')} for \eqn{t=1,\ldots,N}{t=1,...,N},
where \eqn{C} is a scaling constant used
to force the normalization \eqn{\sum{h_{t,k}^2}=1}{sum[t=0,...,N-1]{h(t,k)^2}=1};
\eqn{\widetilde{W}=\pi W (N-1)\Delta t}{W'=pi*W*(N-1)*dt} where
\eqn{\Delta t}{dt} is the sampling interval;
\eqn{g_t=(2t - 1)/N}{g(t)=(2*t-1)/N}; and
\eqn{I_0(\cdot)}{I0()} is the modified Bessel function
of the first kind and zeroth order. The parameter \eqn{W}
is related to the \emph{resolution bandwidth} since
it roughly defines the desired half-width of the central
lobe of the resulting spectral window. Higher order
DPSS tapers (i.e., \eqn{h_{t,k}}{h(t,k)} for \eqn{k > 0}{k > 0})
can be calculated using a relatively simple tridiagonalization
formulation (see the references for more information). Finally,
we note that the sampling interval \eqn{\Delta t}{dt} can
be set to unity without any loss of generality.}}
}
\references{
A. T. Walden, ``Accurate Approximation of a 0th Order Discrete Prolate
Spheroidal Sequence for Filtering and Data Tapering", \emph{Signal Processing},
\bold{18}, 341--8 (1989).

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

S.K.Mitra, J. Kaiser  (1993),
\emph{Handbook for Digital Signal Processing}, John Wiley and Sons, Inc.
}

\seealso{
\code{\link{taper}}.}

\examples{
require(ifultools)
## change plot layout 
gap <- 0.11
old.plt <- splitplot(4,4,1,gap=gap)

## create a plot of all supported tapers and 
## windows 
nms <- c("rectangle", "triangle", "raised cosine",
    "hanning", "hamming", "blackman",
    "nuttall", "gaussian", "kaiser",
    "chebyshev", "born jordan", "sine",
    "parzen", "papoulis", "daniell", "dpss")

for (i in seq(along=nms)){
    if (i > 1) splitplot(4,4,i,gap=gap)
    plot(taper(type=nms[i]))
}

## restore plot layout to initial state 
par(old.plt)
}
\keyword{ts}

