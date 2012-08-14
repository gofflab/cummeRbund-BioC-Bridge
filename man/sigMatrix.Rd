\name{sigMatrix}
\alias{sigMatrix}
\alias{sigMatrix,CuffSet-method}


\title{
sigMatrix
}
\description{
Returns a ggplot2 plot object representing a matrix of significant features. This is a useful synopsis of all significant pairwise comparisons within the dataset.
}
\usage{
\S4method{sigMatrix}{CuffSet}(object,alpha=0.05,level='genes',orderByDist=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class CuffSet.
}
  \item{alpha}{
An alpha value by which to filter multiple-testing corrected q-values to determine significance
}
  \item{level}{
Feature level to be queried for significance (must be one of c('genes','isoforms','TSS','CDS')
}
  \item{orderByDist}{
Logical. If TRUE then samples are re-ordered based on JS-distance from one another (fairly useless unless you have a specific need for this).
}
}
\details{
Creates a matrix plot to illustrate the number of significant features of type 'level' at a given alpha from a cuffdiff run.
}
\value{
A ggplot2 plot object
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Create CuffSet object from sample data
	d<-sigMatrix(a) #Create csDensity plot
	d #Render plot
}
