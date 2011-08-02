\name{csScatter}
\alias{csScatter}

\title{
Scatter Plot
}
\description{
A scatter plot comparing the FPKM  values from two samples in a cuffdiff run.
}
\usage{
csScatter(object, x, y, logMode=TRUE, pseudocount=0.0001, labels, smooth=FALSE, ...)
}

\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet')
}
  \item{x}{
Sample name for x axis
}
  \item{y}{
Sample name for y axis
}
  \item{logMode}{
Logical argument to log2-transform data (default: T )
}
  \item{pseudocount}{
Value to add to zero FPKM values for log transformation (default: 0.0001)
}
  \item{smooth}{
Logical argument to add a smooth-fit regression line
}
  \item{labels}{
A list of tracking_ids or gene_short_names that will be 'callout' points in the plot for reference. Useful for finding genes of interest in the field. Not implemented yet.
}
  \item{\dots}{
Additional arguments to csScatter
}
}
\details{

}
\value{
ggplot object with geom_point and geom_rug layers
}
\references{

}
\author{
Loyal A. Goff
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	a<-readCufflinks() #Create CuffSet object
	genes<-a@genes #Create CuffData object for all genes
	s<-csScatter(genes,'hESC','Fibroblasts',smooth=TRUE) #Create plot object
	s #render plot object
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.