\name{cuffSet-class}
\Rdversion{1.1}
\docType{class}
\alias{cuffSet-class}
\alias{conf_hi,cuffSet-method}
\alias{conf_lo,cuffSet-method}
\alias{csApply,cuffSet-method}
\alias{csBoxplot,cuffSet-method}
\alias{csDensity,cuffSet-method}
\alias{csHist,cuffSet-method}
\alias{csmelt,cuffSet-method}
\alias{expressionBarplot,cuffSet-method}
\alias{expressionPlot,cuffSet-method}
\alias{fpkm,cuffSet-method}
\alias{geneNames,cuffSet-method}
\alias{write.fpkm,cuffSet-method}

\title{Class "cuffSet"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("cuffSet", assayData, phenoData, featureData, fpkm, conf_lo, conf_hi, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{conf_hi}{\code{signature(object = "cuffSet")}: ... }
    \item{conf_lo}{\code{signature(object = "cuffSet")}: ... }
    \item{csApply}{\code{signature(object = "cuffSet")}: ... }
    \item{csBoxplot}{\code{signature(object = "cuffSet")}: ... }
    \item{csDensity}{\code{signature(object = "cuffSet")}: ... }
    \item{csHist}{\code{signature(object = "cuffSet")}: ... }
    \item{csmelt}{\code{signature(object = "cuffSet")}: ... }
    \item{expressionBarplot}{\code{signature(object = "cuffSet")}: ... }
    \item{expressionPlot}{\code{signature(object = "cuffSet")}: ... }
    \item{fpkm}{\code{signature(object = "cuffSet")}: ... }
    \item{geneNames}{\code{signature(object = "cuffSet")}: ... }
    \item{write.fpkm}{\code{signature(object = "cuffSet")}: ... }
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
showClass("cuffSet")
}
\keyword{classes}
