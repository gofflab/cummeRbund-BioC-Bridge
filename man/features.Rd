\name{features}
\alias{features}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Features
}
\description{
Returns a data frame of features from a CuffData, CuffFeatureSet, or CuffFeature object
}
\usage{
features(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet', or 'CuffFeature')
}
}
\details{
If features have been added to cuffData via addFeatures() then these will be presented as well.
}
\value{
A data.frame of feature-level information
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
None
}
\examples{
	\dontrun{
	a<-readCufflinks() #Create CuffSet object
	feature.df<-features(a@genes)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cummeRbund }
\keyword{ accessor }% __ONLY ONE__ keyword per line
