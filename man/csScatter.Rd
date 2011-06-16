\name{csScatter}
\alias{csScatter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Scatter Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
csScatter(object, x, y, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class ('CuffData','CuffFeatureSet')
}
  \item{x}{
Sample name for x axis
}
  \item{y}{
Sample name for y axis
}
  \item{logMode}{
Logical argument to log2-transform data (default: T )
}
  \item{pseudocount}{
Value to add to zero FPKM values for log transformation (default: 0.0001)
}
  \item{smooth}{
Logical argument to add a smooth-fit regression line
}
  \item{labels}{
Not implemented yet
}
  \item{\dots}{
Additional arguments to csScatter
}
}
\details{

}
\value{
ggplot object with geom_point and geom_rug layers
}
\references{

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