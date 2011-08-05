\name{readCufflinks}
\alias{readCufflinks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readCufflinks
}
\description{
This initializes the backend SQLite table and provides a DB connection for all downstream data analysis.
}
\usage{
readCufflinks(dir = getwd(), dbFile = "cuffData.db", geneFPKM = "genes.fpkm_tracking", geneDiff = "gene_exp.diff", isoformFPKM = "isoforms.fpkm_tracking", isoformDiff = "isoform_exp.diff", TSSFPKM = "tss_groups.fpkm_tracking", TSSDiff = "tss_group_exp.diff", CDSFPKM = "cds.fpkm_tracking", CDSExpDiff = "cds_exp.diff", CDSDiff = "cds.diff", promoterFile = "promoters.diff", splicingFile = "splicing.diff", driver = "SQLite", rebuild = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
%%     ~~Describe \code{dir} here~~
}
  \item{dbFile}{
%%     ~~Describe \code{dbFile} here~~
}
  \item{geneFPKM}{
%%     ~~Describe \code{geneFPKM} here~~
}
  \item{geneDiff}{
%%     ~~Describe \code{geneDiff} here~~
}
  \item{isoformFPKM}{
%%     ~~Describe \code{isoformFPKM} here~~
}
  \item{isoformDiff}{
%%     ~~Describe \code{isoformDiff} here~~
}
  \item{TSSFPKM}{
%%     ~~Describe \code{TSSFPKM} here~~
}
  \item{TSSDiff}{
%%     ~~Describe \code{TSSDiff} here~~
}
  \item{CDSFPKM}{
%%     ~~Describe \code{CDSFPKM} here~~
}
  \item{CDSExpDiff}{
%%     ~~Describe \code{CDSExpDiff} here~~
}
  \item{CDSDiff}{
%%     ~~Describe \code{CDSDiff} here~~
}
  \item{promoterFile}{
%%     ~~Describe \code{promoterFile} here~~
}
  \item{splicingFile}{
%%     ~~Describe \code{splicingFile} here~~
}
  \item{driver}{
%%     ~~Describe \code{driver} here~~
}
  \item{rebuild}{
%%     ~~Describe \code{rebuild} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
