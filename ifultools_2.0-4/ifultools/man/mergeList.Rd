%%  WARNING: This file was automatically generated from the associated 
%%  it_util.mid file. Do NOT edit this Rd file to make a change. Instead, 
%%  edit the it_util.mid file in the project MID directory. Once the 
%%  it_util.mid file has been updated, this Rd file, and all other 
%%  documentation (such as corresponding LaTeX, SGML and HTML documentation) 
%%  should be regenerated using the mid.pl Perl script. 

%% R documentation for the mergeList function

\name{mergeList}
\alias{mergeList}
\title{Merges the information from two lists}
\concept{utilities}
\usage{mergeList(x, y)}
\description{Given lists X and Y, this functions replaces
the commonly named objects in X with
those of Y and appends the uncommon Y components to X.
List X is returned as the merged list.}

\arguments{
\item{x}{a list of named objects.}

\item{y}{a list of named objects.}

}
\seealso{
\code{\link{prettyPrintList}}.}

\examples{
## develop lists 
dinner <- list(Entree="Spaghetti and Meatballs",
    Starter="Caesar Salad", Dessert="Spumoni",
    Beverage="Wine and Water")

## oops, we are all out of spumoni and we just got 
## some tiramisu in from the local bakery 
change <- list(Dessert="Tiramisu",Note="Please tip your waiter")

## merge the lists and prett-print 
prettyPrintList(mergeList(dinner, change), header="What's for dinner?")
}
\keyword{list}

