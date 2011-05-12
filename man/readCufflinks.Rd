\name{readCufflinks}
\alias{readCufflinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readCufflinks
}
\description{
Reads cufflinks output file(s) and stores expression information and associated annotation in a cuffSet object.
User can supply additional phenoData and featureData files to extend the default annotation.
}
\usage{
readCufflinks(fpkmFile, 
			phenoDataFile, 
			featureDataFile, 
			experimentDataFile, 
			notesFile, 
			path, 
			annotation, 
			fpkmArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, na.string = na.string, ...), 
			phenoDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, na.string = na.string, stringsAsFactors = stringsAsFactors, ...), 
			featureDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, na.string = na.string, stringsAsFactors = stringsAsFactors, ...), 
			experimentDataArgs = list(sep = sep, header = header, row.names = row.names, quote = quote, na.string = na.string, stringsAsFactors = stringsAsFactors, ...), 
			sep = "\t", 
			na.string = "-", 
			header = TRUE, 
			quote = "", 
			stringsAsFactors = FALSE, 
			row.names = 1L, 
			...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fpkmFile}{
The name of an '*.fpkm_tracking' file output from a cufflinks analysis.
}
  \item{phenoDataFile}{
An optional file with phenotype/parameter information for samples (conditions). (File will have header with parameter names.  First column must match sample names from fpkmFile.)
}
  \item{featureDataFile}{
An optional file with annotation/paramter information for features (genes). (File will have header with parameter names. First column must match feature names from fpkmFile.)
}
  \item{experimentDataFile}{
An optional file with information about experimental conditions. (Not sure about format yet).
}
  \item{notesFile}{
An optional file with misc notes that will be appended to cuffSet object.
}
  \item{path}{
Path to fpkmFile (Defaults to current working directory).
}
  \item{annotation}{
%%     ~~Describe \code{annotation} here~~
}
  \item{fpkmArgs}{
A list of default 'read.table' args for reading fpkmFile (e.g. list(sep = sep, header = header, row.names = row.names, quote = quote, na.string = na.string, ...))
}
  \item{phenoDataArgs}{
A list of 'read.table' args for reading phenoDataFile (see fpkmArgs)
}
  \item{featureDataArgs}{
A list of 'read.table' args for reading featureDataFile (see fpkmArgs)
}
  \item{experimentDataArgs}{
A list of 'read.table' args for reading experimentDataFile (see fpkmArgs)
}
  \item{sep}{
Field separator for input files.
}
  \item{na.string}{
Missing value string for input files.
}
  \item{header}{
T/F input files have header. (Default: T)
}
  \item{quote}{
T/F values are in quotations
}
  \item{stringsAsFactors}{
T/F Treat strings in parameter files as factors (Default: T)
}
  \item{row.names}{
Which row contains row names.
}
  \item{\dots}{
Pass any additional arguments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class cuffSet.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Loyal A. Goff
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cufflinks }
\keyword{ cuffSet }
\keyword{ RNAseq }
\keyword{ nextgen }
