\name{fpkmMatrix}
\alias{fpkmMatrix}
\alias{fpkmMatrix,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fpkmMatrix
}
\description{
Retrieve FPKM values as gene by condition matrix
}
\usage{
\S4method{fpkmMatrix}{CuffData}(object)
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
A feature x condition matrix of FPKM values.
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
	fpkmMatrix(sampleGeneSet)
}
