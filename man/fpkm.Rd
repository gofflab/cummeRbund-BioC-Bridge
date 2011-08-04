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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data.frame of FPKM-level values for a set of features.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
