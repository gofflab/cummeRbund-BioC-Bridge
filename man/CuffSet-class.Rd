\name{CuffSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffSet-class}
\alias{[,CuffSet-method}
\alias{addFeatures,CuffSet-method}
\alias{getGene,CuffSet-method}
\alias{getGenes,CuffSet-method}
\alias{getLevels,CuffSet-method}
\alias{samples,CuffSet-method}

\title{Class "CuffSet"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffSet", DB, conditions, genes, isoforms, TSS, CDS, promoters, splicing, relCDS, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{DB}:}{Object of class \code{"SQLiteConnection"} ~~ }
    \item{\code{conditions}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{genes}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffData"} ~~ }
    \item{\code{promoters}:}{Object of class \code{"CuffDist"} ~~ }
    \item{\code{splicing}:}{Object of class \code{"CuffDist"} ~~ }
    \item{\code{relCDS}:}{Object of class \code{"CuffDist"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "CuffSet")}: ... }
    \item{addFeatures}{\code{signature(object = "CuffSet")}: ... }
    \item{getGene}{\code{signature(object = "CuffSet")}: ... }
    \item{getGenes}{\code{signature(object = "CuffSet")}: ... }
    \item{getLevels}{\code{signature(object = "CuffSet")}: ... }
    \item{samples}{\code{signature(object = "CuffSet")}: ... }
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
showClass("CuffSet")
}
\keyword{classes}
