\name{shannon.entropy}
\alias{shannon.entropy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Shannon entropy
}
\description{
Calculates the Shannon entropy for a probability distribution
}
\usage{
shannon.entropy(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A vector of probabilities (must sum to ~1)
}
}
\details{
None
}
\value{
Returns a numeric value for the Shannon entropy of the supplied probability distribution
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Loyal A. Goff
}
\note{
None
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
	x<-sample(1:500,50)
	p<-x/sum(x)
	shannon.entropy(p)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cummeRbund }
\keyword{ Shannon }% __ONLY ONE__ keyword per line
