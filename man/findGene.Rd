\name{findGene}
\alias{findGene}
\alias{findGene,CuffSet-method}

\title{
findGene
}
\description{
A helper function to retrieve the gene_ids given a 'lookup' value (e.g. gene_short_name, isoform_id, etc). 
Utility to search for gene_id and gene_short_name given a single 'query' string (e.g. query='pink1' will return all genes with 'pink1' (case-insensitive) in the gene_short_name field.
}
\usage{
\S4method{findGene}{CuffSet}(object, query)
}

\arguments{
  \item{object}{
An object of class 'CuffSet' (Primary 'pointer' object for Cufflinks data).
}
  \item{query}{
A character string for which you would like to retrieve corresponding gene_id values.
}
}
\details{
None.
}
\value{
Returns a data.frame of gene_id and gene_short_name values corresponding to genes from which 'query' matches
}
\references{
None.
}
\author{
Loyal A. Goff
}
\note{
Right now, this does not return an error if it cannot find a gene.  (this is probably a bad thing...)
}

\examples{
	cuff<-readCufflinks(system.file("extdata", package="cummeRbund")) #Read cufflinks data and create master CuffSet object
	myQuery<-'pink1'
	findGene(cuff,myQuery) # Retrieve gene_id values for any genes matching 'pink1'
}
