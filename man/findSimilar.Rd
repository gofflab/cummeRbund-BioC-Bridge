\name{findSimilar}
\alias{findSimilar}
\alias{findSimilar,CuffSet-method}

\title{
findSimilar
}
\description{
Returns a CuffGeneSet containing n genes with the most similar expression profiles to gene/profile x.
}
\usage{
\S4method{findSimilar}{CuffSet}(object, x, n)
}

\arguments{
  \item{object}{
A object of class 'CuffSet'
}
  \item{x}{
A 'gene_id' or 'gene_short_name' from which to look up an expression profile OR a vector of expression values to compare all genes (vector must have same length and order of 'samples')
}
  \item{n}{
Number of similar genes to return
}
}
\details{
By default, returns a CuffGeneSet object with n similar genes.  This may change in the future.
}
\value{
A CuffGeneSet object of n most similar genes to x.
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
\examples{
	a<-readCufflinks(system.file("extdata", package="cummeRbund"))
	mySimilarGenes<-findSimilar(a,"PINK1",10)
}