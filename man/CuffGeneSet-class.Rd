\name{CuffGeneSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffGeneSet-class}

\title{Class "CuffGeneSet"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffGeneSet", annotation, fpkm, diff, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{ids}:}{Object of class \code{"character"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffFeatureSet"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CuffFeatureSet}"}, directly.
}
\section{Methods}{
No methods defined with class "CuffGeneSet" in the signature.
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
showClass("CuffGeneSet")
}
\keyword{classes}
