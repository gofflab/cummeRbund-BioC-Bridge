\name{getGene}
\alias{getGene}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getGene
}
\description{
Primary accessor from a CuffSet object to retrive all related information for 1 (one) given gene, indexed by gene_id or gene_short_name.
}
\usage{
getGene(object, geneId, sampleIdList=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class 'CuffSet' (Primary 'pointer' object for Cufflinks data).
}
  \item{geneId}{
A character string to identify which gene for which you would like to retrieve all information.
}
  \item{sampleIdList}{
A vector of sample names used to subset or re-order samples in returned object
}
}
\details{
None.
}
\value{
Returns a CuffGene object containing all related information for a given gene_id or gene_short_name
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
Right now, this does not return an error if it cannot find a gene.  (this is probably a bad thing...)
}

\examples{
	\dontrun{
	a<-readCufflinks() #Read cufflinks data and create master CuffSet object
	myGene<-getGene(a,"PINK1") # Retrieve all information for gene "PINK1"
	}
}
