\name{fpkm}
\alias{fpkm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retrieve FPKM values
}
\description{
Returns a data.frame from @FPKM slot
}
\usage{
fpkm(object, features=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet','CuffFeature')
}
  \item{features}{
A logical argument to merge the @features slot and return a data.frame with all attributes attached.
}
}
\details{
Returns a data.frame of FPKM values.
}
\value{
A data.frame of FPKM-level values for a set of features.
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
	fpkm(PINK1)
}


\keyword{ cummeRbund }
\keyword{ accessor }
