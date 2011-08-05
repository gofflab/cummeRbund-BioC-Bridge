\name{CuffData-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffData-class}
\alias{addFeatures,CuffData-method}
\alias{csBoxplot,CuffData-method}
\alias{csDensity,CuffData-method}
\alias{csScatter,CuffData-method}
\alias{csVolcano,CuffData-method}
\alias{diffData,CuffData-method}
\alias{dim,CuffData-method}
\alias{featureNames,CuffData-method}
\alias{features,CuffData-method}
\alias{fpkm,CuffData-method}
\alias{fpkmMatrix,CuffData-method}
\alias{getFeatures,CuffData-method}
\alias{getLevels,CuffData-method}
\alias{samples,CuffData-method}

\title{Class "CuffData"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffData", DB, tables, filters, type, idField, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{tables}:}{Object of class \code{"list"} ~~ }
    \item{\code{filters}:}{Object of class \code{"list"} ~~ }
    \item{\code{type}:}{Object of class \code{"character"} ~~ }
    \item{\code{idField}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{addFeatures}{\code{signature(object = "CuffData")}: ... }
    \item{csBoxplot}{\code{signature(object = "CuffData")}: ... }
    \item{csDensity}{\code{signature(object = "CuffData")}: ... }
    \item{csScatter}{\code{signature(object = "CuffData")}: ... }
    \item{csVolcano}{\code{signature(object = "CuffData")}: ... }
    \item{diffData}{\code{signature(object = "CuffData")}: ... }
    \item{dim}{\code{signature(x = "CuffData")}: ... }
    \item{featureNames}{\code{signature(object = "CuffData")}: ... }
    \item{features}{\code{signature(object = "CuffData")}: ... }
    \item{fpkm}{\code{signature(object = "CuffData")}: ... }
    \item{fpkmMatrix}{\code{signature(object = "CuffData")}: ... }
    \item{getFeatures}{\code{signature(object = "CuffData")}: ... }
    \item{getLevels}{\code{signature(object = "CuffData")}: ... }
    \item{samples}{\code{signature(object = "CuffData")}: ... }
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
showClass("CuffData")
}
\keyword{classes}
