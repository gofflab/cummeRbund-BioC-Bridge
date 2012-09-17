\name{Dimensionality Reduction}
\alias{MDSplot}
\alias{MDSplot,CuffData-method}
\alias{PCAplot}
\alias{PCAplot,CuffData-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Dimensionality reduction utilities
}
\description{
Dimensionality reduction plots for feature selection and extraction for cummeRbund
}
\usage{
\S4method{MDSplot}{CuffData}(object,replicates=FALSE,logMode=TRUE,pseudocount=1.0)
\S4method{PCAplot}{CuffData}(object,x="PC1", y="PC2",replicates=FALSE,pseudocount=1.0,scale=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The output of class CuffData from which to draw expression estimates. (e.g. genes(cuff))
}
  \item{x}{
For PCAplot, indicates which principal component is to be presented on the x-axis (e.g. "PC1","PC2","PC3", etc) 
}
  \item{y}{
See x.
}
  \item{pseudocount}{
Value added to FPKM to avoid log transformation issues.
}
  \item{logMode}{
Logical value whether or not to use log-transformed expression estimates (default: TRUE)
}
  \item{replicates}{
A logical value to indicate whether or not individual replicate expression estimates will be used.
}
  \item{scale}{
For PCAplot, a logical value passed directly to prcomp.
}
  \item{\dots}{
Additional passthrough arguments (may not be fully implemented yet). 
}
}
\details{
These methods attempt to project a matrix of expression estimates across conditions and/or replicates onto a smaller number of dimesions for feature selection, feature extraction, and can
also be useful for outlier detection.

}		
\value{
A ggplot2 object.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}


\examples{
	cuff<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object from sample data
	p<-PCAplot(genes(cuff),x="PC2",y="PC3",replicates=TRUE)
	m<-MDSplot(genes(cuff),replicates=TRUE)
	p #Render PCA plot
	m #Render MDS plot
	
}

