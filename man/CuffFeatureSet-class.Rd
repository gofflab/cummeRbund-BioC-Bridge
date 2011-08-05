\name{CuffFeatureSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffFeatureSet-class}
\alias{csHeatmap,CuffFeatureSet-method}
\alias{csScatter,CuffFeatureSet-method}
\alias{csVolcano,CuffFeatureSet-method}
\alias{diffData,CuffFeatureSet-method}
\alias{expressionBarplot,CuffFeatureSet-method}
\alias{featureNames,CuffFeatureSet-method}
\alias{features,CuffFeatureSet-method}
\alias{fpkm,CuffFeatureSet-method}
\alias{fpkmMatrix,CuffFeatureSet-method}
\alias{samples,CuffFeatureSet-method}

\title{Class "CuffFeatureSet"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffFeatureSet", annotation, fpkm, diff, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{csHeatmap}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{csScatter}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{csVolcano}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{diffData}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{expressionBarplot}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{featureNames}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{features}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{fpkm}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{fpkmMatrix}{\code{signature(object = "CuffFeatureSet")}: ... }
    \item{samples}{\code{signature(object = "CuffFeatureSet")}: ... }
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
showClass("CuffFeatureSet")
}
\keyword{classes}
