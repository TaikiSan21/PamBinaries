%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the isVectorAtomic function

\name{isVectorAtomic}
\alias{isVectorAtomic}
\title{Tests whether an object is a vector or not}
\concept{utilities}
\usage{isVectorAtomic(x)}
\description{The \code{is.vector} function returns a \code{FALSE}
value in some cases where intuitively one might
expect a \code{TRUE} value to be returned. For example,
\code{is.vector(z)} returns \code{FALSE} for each of the
following:
\describe{
\item{i}{\code{z <- 1:3;names(z) <- 1:3}}

\item{ii}{\code{z <- matrix(1:3, nrow=1)}}

\item{iii}{\code{z <- matrix(1:3, ncol=1)}}}

These results are not necessarily incorrect, they are
just one interpretion of the definition of a vector.
Contrarily, the \code{isVectorAtomic(z)} function returns
\code{TRUE} for each of the above examples. Thus,
\code{isVectorAtomic} expands the basic definition of a
vector to allow matrices containing a single row
or column and named vectors. Also, unlike \code{is.vector},
\code{isVectorAtomic} returns \code{FALSE} for objects of class \code{list}.}

\arguments{
\item{x}{an object of arbitrary class.}

}

\value{
a vector of character strings containing the result. The length
of this vector is equal to \code{length(x)}.
}
\seealso{
\code{\link{rotateVector}}, \code{\link{mergeList}}.}

\examples{
## cases where isVectorAtomic returns TRUE 
z <- 1:3;names(z) <- letters[1:3]
isVectorAtomic(z)
isVectorAtomic(matrix(1:3, nrow=1))
isVectorAtomic(matrix(1:3, ncol=1))
isVectorAtomic(1:5)
isVectorAtomic(letters)

## cases where isVectorAtomic returns FALSE 
isVectorAtomic(list(1:3))
isVectorAtomic(data.frame(1:3,2:4))
isVectorAtomic(matrix(1:12, nrow=4))
}
\keyword{utilities}

