\name{getSig}
\alias{getSig}
\alias{getSig,CuffSet-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getSig
}
\description{
Returns the identifiers of significant genes in either a list format (by pairwise comparison) or in a test table - like format.
}
\usage{
getSig(object,x,y,level='genes',testTable=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A CuffSet object (e.g. cuff)
}
  \item{x}{
Not implemented yet.
}
  \item{y}{
See x.
}
  \item{level}{
Feature level to be queried for significance (must be one of c('genes','isoforms','TSS','CDS')
}
  \item{testTable}{
Logical value whether to return the results as a list of ID sets (default), or a test table of binary significance results (1=significant, 0=not significant).
}
}
\details{
This is a convenience function to quickly retrieve lists of identifiers for genes or features that were determined to be significantly regulated between conditions by cuffdiff.
This function relies (for now) on the 'significant' field of the cuffdiff diffData files. Future versions may allow for filtering significance lists by q-value instead.
By default getSig returns a list of ID-sets, where the names for the list represent the pairwise comparisons from the cuffdiff analysis.  If testTable=T then getSig will ouput a
'test table' of significant genes X pairwise tests, where the values in the table are 1 for features that are significant for this comparison and 0 for genes that are not. TestTable=T 
only includes the features that are significant in at least one comparison.
}		
\value{
Either a list of feature ID sets, or a data.frame of pairwise test results.
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
None.
}


\examples{

	a<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data in sample directory and creates CuffSet object 'a'
	mySig<-getSig(a,level='genes')
	myTestTable<-getSig(a,level='genes',testTable=TRUE)
	names(mySig)
	head(myTestTable)
	
}

