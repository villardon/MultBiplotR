\name{GowerProximities}
\alias{GowerProximities}
\title{Gower Dissimilarities for mixed types of data}
\description{
  Gower Dissimilarities for mixed types of data
}
\usage{
  GowerProximities(x, y = NULL, Binary = NULL, Classes = NULL,
                   transformation = 3, IntegerAsOrdinal = FALSE, BinCoef
                   = "Simple_Matching", ContCoef = "Gower", NomCoef =
                   "GOW", OrdCoef = "GOW")
}

\arguments{
 \item{x}{
Main data. Distances among rows are calculated if y=NULL. Must be a data frame.
}
  \item{y}{
Suplementary data matrix. If not NULL the distances among the rows of x and y are calculated. Must be a data frame with the same columns as x.
}
  \item{Binary}{
A vector containing the binary variables.
}
  \item{Classes}{
Vector with column types. If NULL the data frame types are used.
}
  \item{transformation}{
Transformation for the similarities.
}
  \item{IntegerAsOrdinal}{
Should integer variables be used as ordinal?
}
  \item{BinCoef}{
Coefficient for the binary data
}
  \item{ContCoef}{
Coefficient for the continuous data
}
  \item{NomCoef}{
Coefficient for the nominal data
}
  \item{OrdCoef}{
Coefficient for the ordinal data
}
}

\details{

The transformation \code{sqrt(1-S)} is applied to the similarity.
}
\value{
An object of class \code{proximities}.This has components:
\item{comp1 }{Description of }

}
\references{
J. C. Gower. (1971) A General Coefficient of Similarity and Some of its Properties. Biometrics, Vol. 27, No. 4, pp. 857-871.
}
\author{
Jose Luis Vicente-Villardon
}
\examples{
data(spiders)
}
