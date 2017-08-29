%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the prettyPrintList function

\name{prettyPrintList}
\alias{prettyPrintList}
\title{Prints a list of named objects}
\concept{utilities}
\usage{prettyPrintList(x, header=NULL, justify="left", sep=":")}
\description{Prints a list of named objects in a two column format, with the
names of objects in the list in the first column and a summary of the object
in the second. The user can control the justification and separation character
between the columns.}

\arguments{
\item{x}{a list of named objects.}
\item{header}{a character string defining the header. Default: \code{NULL} (no header).}

\item{justify}{a character string defining the aligment format for the
objects in the list. Choices are \code{"none"}, \code{"left"}, \code{"right"} and \code{"decimal"}.}

\item{sep}{a siingle character used to divide the list object names from their
corresponding values.}
}
\seealso{
\code{\link{isVectorAtomic}}.}

\examples{
## develop menu 
dinner <- list(Entree="Spaghetti and Meatballs",
    Starter="Caesar Salad", Dessert="Spumoni",
    Beverage="Wine and Water")

## pretty-print the menu 
prettyPrintList(dinner, header="What's for dinner?")
}
\keyword{list}

