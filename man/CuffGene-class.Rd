\name{CuffGene-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffGene-class}

\title{Class "CuffGene"}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffGene", annotation, fpkm, diff, ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} ~~ }
    \item{\code{isoforms}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{TSS}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{CDS}:}{Object of class \code{"CuffFeature"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{fpkm}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{diff}:}{Object of class \code{"data.frame"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{CuffFeature}"}, directly.
}
\section{Methods}{
No methods defined with class "CuffGene" in the signature.
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
showClass("CuffGene")
}
\keyword{classes}
