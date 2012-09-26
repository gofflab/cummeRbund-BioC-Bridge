\name{features}
\alias{features}
\alias{features,CuffGene-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Features
}
\description{
Returns a data frame of features from a CuffGene object
}
\usage{
\S4method{features}{CuffGene}(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffGene')
}
}
\details{
None
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

\examples{
	data(sampleData)
	features(PINK1)
}
