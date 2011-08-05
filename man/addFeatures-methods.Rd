\name{addFeatures-methods}
\docType{methods}
\alias{addFeatures-methods}
\alias{addFeatures,CuffData-method}
\alias{addFeatures,CuffSet-method}
\title{ ~~ Methods for Function addFeatures in Package `cummeRbund' ~~}
\description{
 ~~ Methods for function \code{addFeatures} in Package `cummeRbund' ~~
}
\section{Methods}{
\describe{

\item{\code{signature(object = "CuffData")}}{
Adds a data.frame of feature-level data (first column must be 'id' value for that feature-type [e.g. 'gene_id' for genes, 'isoform_id' for isoforms, etc]) to cuffData.db database
}

\item{\code{signature(object = "CuffSet")}}{
Adds a data.frame of feature-level data (first column must be 'id' value for that feature-type [e.g. 'gene_id' for genes, 'isoform_id' for isoforms, etc]) to cuffData.db database. Defaults to gene-level features unless otherwise specified.
}
}}
\keyword{methods}
\keyword{setter}
