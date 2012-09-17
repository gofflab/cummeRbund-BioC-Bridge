\name{QCplots}
\alias{fpkmSCVPlot}
\alias{fpkmSCVPlot,CuffData-method}

\title{
Quality Control visualizations		
}
\description{
A collection of ggplot2 visualizations for quality control assessment of cuffdiff output.

- fpkmSCVPlot: A measure of cross-replicate variability, the squared coefficient of variation is a normalized measure of variance between empirical repicate FPKM values per condition, across the range of FPKM estimates.

}
\usage{
\S4method{fpkmSCVPlot}{CuffData}(object,FPKMLowerBound=1)
}

\arguments{
  \item{object}{
An object of class CuffData.
}
  \item{FPKMLowerBound}{
A lower limit cutoff for FPKM values from which a fit of squared Coefficient of variation (default: 1)  
}
}
\details{
None
}
\value{
A ggplot2 plot object with a geom_box layer.
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
	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create CuffSet object
	genes<-a@genes #CuffData object for all genes
	csBoxplot(genes)
}

