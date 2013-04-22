\name{Exploratory Analysis}
\alias{csNMF}
\alias{csNMF,CuffData-method}
\alias{csNMF,CuffFeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for Exploratory Analysis in cummeRbund
}
\description{
Exploratory analysis methods for cummeRbund RNA-Seq data.
}
\usage{
\S4method{csNMF}{CuffData}(object,k,logMode=T,pseudocount=1,maxiter=1000,replicates=FALSE,fullnames=FALSE)
\S4method{csNMF}{CuffFeatureSet}(object,k,logMode=T,pseudocount=1,maxiter=1000,replicates=FALSE,fullnames=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
The output of class CuffData or CuffFeatureSet from which to draw expression estimates. (e.g. genes(cuff) or custom feature set via getGenes() or getFeatures() )
}
  \item{k}{
rank value for factorization
}
  \item{logMode}{
Logical value whether or not to use log-transformed FPKM values. [Default: TRUE]
}
  \item{pseudocount}{
Value added to FPKM to avoid log transformation issues.
}
  \item{maxiter}{
Maximum number of iterations for factorization [Default: 1000]
}
  \item{replicates}{
A logical value to indicate whether or not individual replicate expression estimates will be used.
}
  \item{fullnames}{
Logical passthrough value to fpkmMatrix whether or not to concatenate gene_short_name with tracking_id. [Default: FALSE]
}
}
\details{
csNMF is a convenience method to invoke the nnmf() method from package:NMFN.  This performs non-negative matrix factorization on the provided data and can be useful for many downstream applications.

}		
\value{
csNMF returns W, H - decomposed matrices of input FPKM values. (See package:NMFN for details)
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
	data(sampleData)
	csNMF(sampleGeneSet,4)
	
}

