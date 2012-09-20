\name{csDistHeat}
\alias{csDistHeat}
\alias{csDistHeat,CuffFeatureSet-method}
\alias{csDistHeat,CuffData-method}

\title{
  csDistHeat
}
\description{
  Creates a ggplot plot object with a geom_tile layer of JS Distance
  values between samples or genes.
}
\usage{
\S4method{csDistHeat}{CuffFeatureSet}(object, replicates=F, samples.not.genes=T,
  logMode=T, pseudocount=1.0,
  heatscale=c(low='lightyellow',mid='orange',high='darkred'),
  heatMidpoint=NULL, ...)
}

\arguments{
  \item{object}{
    An object of class 'CuffFeatureSet' or 'CuffGeneSet'
  }
   \item{replicates}{
    A logical argument whether or not to use individual replicate FPKM values as opposed to condition FPKM estimates. (default: FALSE)
  }
  \item{samples.not.genes}{
    Compute distances between samples rather than genes. If False,
    compute distances between genes.
  }
  \item{logMode}{
    A logical argument to log10-transform FPKM values prior to plotting.
  }
  \item{pseudocount}{
    Value to be added to FPKM for appropriate log transformation and clustering. (Avoids zero-based errors)
  }
  \item{heatscale}{
    A list with min length=2, max length=3 that describe the the color scale.
  }
  \item{heatMidpoint}{
    Value for midpoint of color scale.
  }
  \item{\dots}{
    Additional arguments to csHeatmap
  }
}
\details{
  None
}
\value{
  A ggplot2 plot object with a geom_tile layer to display distance
  between samples or genes.
}
\references{
  None
}
\author{
  Loyal A. Goff, Cole Trapnell, and David Kelley
}
\note{
  None
}

\examples{
	data(sampleData)
	csDistHeat(sampleGeneSet) 
}

\keyword{heatmap}
\keyword{Jensen-Shannon}
\keyword{distance}
