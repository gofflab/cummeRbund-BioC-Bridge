\name{csBoxplot}
\alias{csBoxplot}

\title{
csBoxplot
}
\description{
Creates a ggplot2 plot object with a geom_box layer displaying summary statistics for FPKM values across samples (x).
}
\usage{
csBoxplot(object, logMode=T, ...)
}

\arguments{
  \item{object}{
An object of class CuffData.
}
  \item{logMode}{
A logical argument to log10 -transform FPKM values.
}
  \item{\dots}{
Additional arguments to csBoxplot
}
}
\details{
None
}
\value{
A ggplot2 plot object with a geom_box layer.
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
	\dontrun{
	a<-readCufflinks #Read cufflinks data and create CuffSet object
	genes<-a@genes #CuffData object for all genes
	csBoxplot(genes)
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
