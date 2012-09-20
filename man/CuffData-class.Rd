\name{CuffData-class}
\Rdversion{1.1}
\docType{class}
\alias{CuffData-class}
\alias{dim,CuffData-method}
\alias{getFeatures,CuffData-method}
\alias{DB,CuffData-method}
\alias{diffTable}
\alias{makeRnk}
\alias{diffTable,CuffData-method}
\alias{makeRnk,CuffData-method}
\alias{annotation,CuffData-method}

\title{Class "CuffData"}
\description{
A 'pointer' class for all information (FPKM, annotation, differential expression) for a given feature type (genes, isoforms, TSS, CDS). The methods for this function
communicate directly with the SQL backend to present data to the user.
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
    \item{dim}{\code{signature(x = "CuffData")}: ... }
    \item{getFeatures}{\code{signature(object = "CuffData")}: ... }
    \item{DB}{\code{signature(object = "CuffData")}: Accessor for @DB slot }
    \item{diffTable}{\code{signature(object = "CuffData")}: Create a Full table (wide format) of differential expression information for all pairwise comparisons }
    \item{makeRnk}{\code{signature(object = "CuffData")}: Internal method to create .rnk file. Should not be called directly }	
    \item{annotation}{\code{signature(object="CuffData")}: Access annotation data}
	 }
}
\references{
None
}
\author{
Loyal A. Goff
}
\note{
None
}



\seealso{
None
}
\examples{
showClass("CuffData")
}
\keyword{classes}
