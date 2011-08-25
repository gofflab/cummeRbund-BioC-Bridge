\name{createDB}
\alias{createDB}

\title{
createDB
}
\description{
This should not be called directly by the user.
}
\usage{
createDB(dbFname = "cuffData.db", driver = "SQLite")
}

\arguments{
  \item{dbFname}{
File name for backend database (by default this is 'cuffData.db').  If you change this value, be sure to point to the new file with every call to 'readCufflinks'.
}
  \item{driver}{
DB driver for backend database (only SQLite supported at this time.)
}
}
\details{
Should not be called directly
}
\value{
Creates database backend file
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
