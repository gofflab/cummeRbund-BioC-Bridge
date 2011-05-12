\name{JSdist}
\alias{JSdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jensen-Shannon distance on columns
}
\description{
JSdist takes a matrix of expression probabilites (calculated directly or output from makeprobs()) and returns a dist object 
of the pairwise Jensen-Shannon distances between columns
}
\usage{
JSdist(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{
A matrix of expression probabilities (e.g. from makeprobs())
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A dist object of pairwise J-S distances between columns.
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
library(cummeRbund)

mat<-matrix(sample(1:50,50),10)
probs<-makeprobs(mat)
js.distance<-JSdist(probs)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Jensen-Shannon }
\keyword{ distance }
\keyword{ metric }
\keyword{ JS }