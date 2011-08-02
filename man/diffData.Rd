\name{diffData}
\alias{diffData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential comparison data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
diffData(object, x, y, features=FALSE, lnFcCutoff=20, ...)
}

\arguments{
  \item{object}{
An object of class ('CuffData' or 'CuffFeatureSet')
}
  \item{x}{
  Optional, if x and y are both missing, data for all pairwise differential comparisons are returned, otherwise if x and y are sample names from the 'samples' table,
  than only differential data pertaining to those two samples are returned.
}
  \item{y}{
See 'x'
}
  \item{features}{
A logical value that returns all feature-level data as part of data.frame when true.  object must be of class 'CuffData'.
}
  \item{lnFcCutoff}{
For plotting purposes, extreme values of ln_fold_change are removed from view. This is the cutoff for where these should be trimmed.
}

  \item{\dots}{
Additional arguments.
}
}
\details{

}
\value{
A data.frame object
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

}
\examples{
	a<-readCufflinks() #Create CuffSet object
	diff<-diffData(a@genes) #returns a dataframe of differential expression data for all genes.
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cummeRbund }
\keyword{ differential expression }% __ONLY ONE__ keyword per line
