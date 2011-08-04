\name{expressionPlot}
\alias{expressionPlot}

\title{
Expression Plot
}
\description{
A line plot (optionally with confidence intervals) detailing FPKM expression levels across conditions for a given gene(s) or feature(s)
}
\usage{
expressionPlot(object, logMode=FALSE, pseudocount=0.0001, drawSummary=FALSE, sumFun=mean_cl_boot, showErrorbars=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffFeature' or 'CuffGene')
}
  \item{logMode}{
A logical value to draw y-axis (FPKM) on log-10 scale. Default = FALSE.
}
  \item{pseudocount}{
A numeric value added to FPKM to avoid errors on log-10 transformation.
}
  \item{drawSummary}{
A logical valuie. Draws a 'summary' line with mean FPKM values for each condition.
}
  \item{sumFun}{
Function used to determine values for summary line. Default = mean_cl_boot
}
  \item{showErrorbars}{
A logical value whether or not to draw error bars.
}
  \item{\dots}{
Additional arguments
}
}
\details{
None
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



\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	\dontrun{
	a<-readCufflinks() #Create CuffSet object
	myGeneId<-"gene_x" 
	myGene<-getGene(a,myGeneId) #Create CuffFeature object
	expressionPlot(myGene) #Line plot of myGene FPKM values
	expressionPlot(myGene@isoforms) #Line plot of myGene faceted by isoforms
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
