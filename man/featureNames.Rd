\name{featureNames}
\alias{featureNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Feature names
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
featureNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData' or 'CuffFeatureSet')
}
}
\details{
None
}
\value{
A list of feature names
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
	names<-featureNames(a@genes)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cummeRbund }
\keyword{ accessor }% __ONLY ONE__ keyword per line
