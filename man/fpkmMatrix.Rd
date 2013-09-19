\name{fpkmMatrix}
\alias{fpkmMatrix}
\alias{fpkmMatrix,CuffData-method}
\alias{fpkmMatrix,CuffFeatureSet-method}
\alias{repFpkmMatrix}
\alias{repFpkmMatrix,CuffData-method}
\alias{repFpkmMatrix,CuffFeatureSet-method}
\alias{repFpkmMatrix,CuffFeature-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
fpkmMatrix
}
\description{
Retrieve FPKM values as gene by condition (fpkmMatrix) or gene by replicate (repFpkmMatrix) matrix
}
\usage{
\S4method{fpkmMatrix}{CuffData}(object,fullnames=FALSE,sampleIdList)
\S4method{repFpkmMatrix}{CuffData}(object,fullnames=FALSE,repIdList)
}

\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet','CuffGeneSet','CuffGene',or 'CuffFeature')
}
  \item{fullnames}{
A logical value whether or not to concatenate gene_short_name and tracking_id values (easier to read labels)
}
  \item{sampleIdList}{
A vector of sample names to subset the resulting matrix.
}
  \item{repIdList}{
A vector of sample names to subset the resulting matrix.
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
	repFpkmMatrix(sampleGeneSet)
}
