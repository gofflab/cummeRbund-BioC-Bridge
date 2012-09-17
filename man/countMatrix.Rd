\name{countMatrix}
\alias{countMatrix}
\alias{countMatrix,CuffData-method}
\alias{countMatrix,CuffFeatureSet-method}
\alias{repCountMatrix}
\alias{repCountMatrix,CuffData-method}
\alias{repCountMatrix,CuffFeatureSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
countMatrix
}
\description{
Retrieve count values as gene by condition matrix
}
\usage{
\S4method{countMatrix}{CuffData}(object,fullnames=FALSE,sampleIdList)
\S4method{repCountMatrix}{CuffData}(object,fullnames=FALSE,repIdList)
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
A vector of replicate names to subset the resulting replicate matrix.
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
	repCountMatrix(sampleGeneSet)
}
