\name{Convert2ThreeWay}
\alias{Convert2ThreeWay}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three way array from a two way matrix
}
\description{
Converts a two-dimensional matrix into a list where
each cell is the two dimensional data matrix  
for an occasion or group.

}
\usage{
Convert2ThreeWay(x, groups, columns = FALSE, RowNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The two dimensional matrix
}
  \item{groups}{
A factor defining the groups
}
  \item{columns}{
Are the grouos defined for columns?
}
  \item{RowNames}{
Names for the rows of each table.
}
}
\details{
Converts a two dimensional matrix into a multitable list according to the groups provided by the user. Each field of the list has the name of the corresponding group.
}
\value{
A Multitable list. Ech filed is the data matrix for a group.
\item{X }{The multitable list}
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(Chemical)
x= Chemical[,5:16]
X=Convert2ThreeWay(x,Chemical$WEEKS, columns=FALSE)
}