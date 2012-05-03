\name{repFpkmMatrix}
\alias{repFpkmMatrix}
\alias{repFpkmMatrix,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
repFpkmMatrix
}
\description{
Retrieve repFpkm values as gene by condition matrix
}
\usage{
\S4method{repFpkmMatrix}{CuffData}(object)
}

\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet','CuffGeneSet','CuffGene',or 'CuffFeature')
}
}
\details{
None.
}
\value{
A feature x condition matrix of replicate FPKM values.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None
}

\examples{
	data(sampleData)
	repFpkmMatrix(sampleGeneSet)
}
