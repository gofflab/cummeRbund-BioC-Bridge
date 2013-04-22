\name{csCluster}
\alias{csCluster}
\alias{csCluster,CuffFeatureSet-method}

\title{
csCluster
}
\description{
Returns a ggplot2 plot object with geom_line layer plotting FPKM values over conditions faceted by k-means clustering clusters. (Euclidean).  This is very crude at this point.
This does not return any of the clustering information directly, but if you want it, you can retrieve it from the ggplot object returned.
}
\usage{
\S4method{csCluster}{CuffFeatureSet}(object,k,logMode=T,method = "none",pseudocount=1,...)
}
\arguments{
  \item{object}{
An object of class CuffFeatureSet.
}
  \item{k}{
Number of pre-defined clusters to attempt to find.
}
  \item{logMode}{
A logical value whether or not to log-transform the FPKM values prior to clustering.
}
  \item{method}{
Distance function to use when computing cluster solution.  Default "none" will use the Jensen-Shannon distance (JSdist). Provide a function that returns a dist object on rows.
}
  \item{pseudocount}{
Value added to FPKM to avoid log-transform issues.
}
  \item{\dots}{
Additional arguments to pam.
}
}
\details{
Uses 'kmeans' function.
}
\author{
Loyal A. Goff
}
\source{
None
}
\references{
None.
}
\examples{
	data(sampleData)
	csCluster(sampleGeneSet,4)
}
\keyword{datasets}
