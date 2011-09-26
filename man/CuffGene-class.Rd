\name{CuffGene-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffGene-class}
\alias{genes,CuffGene-method}
\alias{isoforms,CuffGene-method}
\alias{TSS,CuffGene-method}
\alias{CDS,CuffGene-method}

\title{Class "CuffGene"}
\description{
A 'data' container class for all FPKM, annotation, and differential expression Data (as well as for all linked features) for a given gene.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CuffGene", annotation, fpkm, diff, ...)}.

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

\section{Accessors}{
	\describe{
		\item{genes}{\code{signature(object = "CuffGene")}: Access @genes slot }
		\item{isoforms}{\code{signature(object = "CuffGene")}: Access @isoforms slot }
		\item{TSS}{\code{signature(object = "CuffGene")}: Access @TSS slot }
		\item{CDS}{\code{signature(object = "CuffGene")}: Access @CDS slot }
	}
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}



\seealso{
\code{\linkS4class{CuffFeature}}
}
\examples{
showClass("CuffGene")
}
\keyword{classes}
