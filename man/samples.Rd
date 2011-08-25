\name{samples}
\alias{samples}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get sample list from CuffData object
}
\description{
Returns a list of sample names from a CuffData or CuffFeatureSet object
}
\usage{
samples(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet','CuffFeature')
}
}
\details{
None
}
\value{
A list of sample names
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
	samples(a@gene)
	}
}

