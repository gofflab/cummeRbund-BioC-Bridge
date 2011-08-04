\name{makeprobs}
\alias{makeprobs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform a matrix into probabilities by columns
}
\description{
This function takes a matrix of expression values (must be greater than 0) and returns a matrix of probabilities by column.
This is a required transformation for the Jensen-Shannon distance which is a metric that operates on probabilities.
}
\usage{
makeprobs(a)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{
A matrix of expression values (values must be greater than 0).
}
}
\details{
To make a matrix of probabilities by row, use t() to transpose prior to calling makeprobs.
}
\value{
A matrix of expression probabilities by column.
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
myMat<-matrix(sample(1:50,50),10)
probs<-makeprobs(myMat)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ probabilities }
\keyword{ transformation }% __ONLY ONE__ keyword per line
