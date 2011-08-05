\name{CuffDist-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffDist-class}
\alias{dim,CuffDist-method}
\alias{distValues,CuffDist-method}

\title{Class "CuffDist"}
\description{
A 'pointer' class to information relative to the distribution-level tests (promoters, splicing, and relative CDS usage)
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffDist", DB, table, type, testId, ...)}.
 
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{table}:}{Object of class \code{"character"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{testId}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "CuffDist")}: ... }
    \item{distValues}{\code{signature(object = "CuffDist")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("CuffDist")
}
\keyword{classes}
