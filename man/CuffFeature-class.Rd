\name{CuffFeature-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffFeature-class}
\alias{expressionBarplot,CuffFeature-method}
\alias{expressionPlot,CuffFeature-method}
\alias{fpkm,CuffFeature-method}
\alias{length,CuffFeature-method}

\title{Class "CuffFeature"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffFeature", annotation, fpkm, diff, ...)}.
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
    \item{expressionBarplot}{\code{signature(object = "CuffFeature")}: ... }
    \item{expressionPlot}{\code{signature(object = "CuffFeature")}: ... }
    \item{fpkm}{\code{signature(object = "CuffFeature")}: ... }
    \item{length}{\code{signature(x = "CuffFeature")}: ... }
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
showClass("CuffFeature")
}
\keyword{classes}
