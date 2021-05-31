\name{CategoricalProximities}
\alias{CategoricalProximities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Proximities among individuals using nominal variables.
}
\description{
Proximities among individuals using nominal variables.
}
\usage{
CategoricalProximities(Data, SUP = NULL, coefficient = "GOW", transformation = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
A data frame containing categorical (nominal) variables
}
  \item{SUP}{
Supplementary data (Used to project supplementary individuals onto the PCoA configuration, for example)
}
  \item{coefficient}{
Similarity coefficient to use (see details)
}
  \item{transformation}{
Transformation of the similarity into a distance
}
  \item{\dots}{
Extra parameters
}
}
\details{
The function calculates similarities and dissimilarities among a set ob ogjects characterized by a set of nominal variables. The function uses similarities and converts into dissimilarities using a variety of transformations controled by the user.
}
\value{
A list of Values
}
\references{
dos Santos, T. R., & Zarate, L. E. (2015). Categorical data clustering: What similarity measure to recommend?. Expert Systems with Applications, 42(3), 1247-1260.

Boriah, S., Chandola, V., & Kumar, V. (2008). Similarity measures for categorical data: A comparative evaluation. red, 30(2), 3.

}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{data(Doctors)
Dis=CategoricalProximities(Doctors, SUP=NULL, coefficient="GOW" , transformation=3)
pco=PrincipalCoordinates(Dis)
plot(pco, RowCex=0.7, RowColors=as.integer(Doctors[[1]]), RowLabels=as.character(Doctors[[1]]))}
}
