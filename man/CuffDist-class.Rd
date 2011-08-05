\name{CuffDist-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffDist-class}
\alias{dim,CuffDist-method}
\alias{values,CuffDist-method}

\title{Class "CuffDist"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffDist", DB, table, type, testId, ...)}.
%%  ~~ describe objects here ~~ 
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
    \item{values}{\code{signature(object = "CuffDist")}: ... }
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
