\name{CategoricalDistances}
\alias{CategoricalDistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distances among individuals using nominal variables.
}
\description{
Distances among individuals using nominal variables.
}
\usage{
CategoricalDistances(x, y = NULL, coefficient = "GOW", transformation = "sqrt(1-S)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Matrix of Categorical Data
}
\item{y}{
A second matrix of categorical data with the same variables as x
}
\item{coefficient}{
Similarity coefficient to use (see details)
}
  \item{transformation}{
Transformation of the similarity into a distance
}

}
\details{
The function calculates similarities and dissimilarities among a set ob ogjects characterized by a set of nominal variables. The function uses similarities and converts into dissimilarities using a variety of transformations controled by the user.
}
\value{
A matrix with distances among the rows of x and y. If y is NULL the interdistances among the rows of x are calculated.
}
\references{
dos Santos, T. R., & Zarate, L. E. (2015). Categorical data clustering: What similarity measure to recommend?. Expert Systems with Applications, 42(3), 1247-1260.

Boriah, S., Chandola, V., & Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation. red, 30(2), 3.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
##---- Should be DIRECTLY executable !! ----
}