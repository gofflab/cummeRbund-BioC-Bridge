\name{countMatrix}
\alias{countMatrix}
\alias{countMatrix,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
countMatrix
}
\description{
Retrieve count values as gene by condition matrix
}
\usage{
\S4method{countMatrix}{CuffData}(object)
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
A feature x condition matrix of count values.
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
	countMatrix(sampleGeneSet)
}
