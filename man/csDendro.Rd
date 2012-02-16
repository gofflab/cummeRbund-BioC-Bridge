\name{csDendro}
\alias{csDendro}
\alias{csDendro,CuffFeatureSet-method}
\alias{csDendro,CuffData-method}

\title{
csDendro
}
\description{
Creates a grid graphics plot of a dendrogram of Jensen-Shannon distances between conditions of a CuffFeatureSet or CuffGeneSet object.
}
\usage{
\S4method{csDendro}{CuffFeatureSet}(object,logMode=T,pseudocount=1)
\S4method{csDendro}{CuffData}(object,logMode=T,pseudocount=1)
}

\arguments{
  \item{object}{
An object of class 'CuffFeatureSet' or 'CuffGeneSet'
}
  \item{logMode}{
A logical argument to log10-transform FPKM values prior to plotting.
}
  \item{pseudocount}{
Value to be added to FPKM for appropriate log transformation and clustering. (Avoids zero-based errors)
}
 
}
\details{
None
}
\value{
Returns a dendrogram object and plots that object by default.
}
\references{
None.
}
\author{
Loyal A. Goff and Cole Trapnell
}
\note{
None
}

\examples{
	data(sampleData)
	csDendro(sampleGeneSet)
}

\keyword{heatmap}
