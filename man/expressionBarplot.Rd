\name{expressionBarplot}
\alias{expressionBarplot}

\title{
Barplot
}
\description{
A barplot of FPKM values with confidence intervals for a given gene, set of genes, or features of a gene (e.g. isoforms, TSS, CDS, etc).
}
\usage{
expressionBarplot(object, logMode=FALSE, pseudocount=0.0001, showErrorbars=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffFeatureSet','CuffGeneSet','CuffFeature','CuffGene')
}
  \item{logMode}{
A logical value whether or not to draw y-axis on log10 scale. Default = FALSE.
}
  \item{pseudocount}{
Numerical value added to each FPKM during log-transformation to avoid errors.
}
  \item{showErrorbars}{
A logical value whether or not to draw error bars. Default = TRUE
}
  \item{\dots}{
Additional arguments.
}
}
\details{

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
Need to implement logMode and features for this plotting method.
}



\seealso{

}
\examples{
	a<-readCufflinks() #Create CuffSet object
	myGeneId<-"gene_x" 
	myGene<-getGene(a,myGeneId) #Create CuffFeature object
	expressionBarplot(myGene) #Barplot of myGene FPKM values
	expressionBarplot(myGene@isoforms) #Barplot of myGene faceted by isoforms
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cummeRbund }
\keyword{ plotting }
\keyword{ barplot }
