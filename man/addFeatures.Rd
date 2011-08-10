\name{addFeatures}
\alias{addFeatures}

\title{
addFeatures
}
\description{
Adds a data.frame of features to a the SQLite backend database.
}
\usage{
addFeatures(object, features, ...)
}

\arguments{
  \item{object}{
An object of class ('CuffSet' or 'CuffData')
}
  \item{features}{
A data.frame of features to add.  1st column MUST contain ids (ie. gene_id for 'gene' features, isoform_id for 'isoform' features, etc)
}
  \item{\dots}{
Additional arguments.
}
}
\details{
None
}
\value{
None
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
	features #A data.frame of features with gene_id as first column 
	a<-readCufflinks() #Read cufflinks data and create CuffSet object
	addFeatures(a,features) #Adds features to table in cuffData.db for use by CuffData objects.
	}
}
