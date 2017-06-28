\name{GowerSimilarities}
\alias{GowerSimilarities}
\title{Gower Dissimilarities for mixed types of data}
\description{
  Gower Dissimilarities for mixed types of data
}
\usage{
  GowerSimilarities(x, y=NULL, transformation="sqrt(1-S)") 
}

\arguments{
 \item{x}{
Main data. Distances among rows are calculated if y=NULL. Must be a data frame.
}
  \item{y}{
Suplementary data matrix. If not NULL the distances among the rows of x and y are calculated. Must be a data frame with the same columns as x.
}
  \item{transformation}{
Transformatio to apply to the similarities.
}
}

\details{
Gower Dissimilarities for mixed types of data.
The  transformation \code{sqrt(1-S)} is applied to the similarity by default.
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
