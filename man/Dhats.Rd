\name{Dhats}
\alias{Dhats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of Disparities
}
\description{
Calculation of Disparities for a MDS model
}
\usage{
Dhats(P, D, W, Model = c("Identity", "Ratio", "Interval", "Ordinal"), Standardize = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
A matrix of proximities (usually dissimilarities)
}
  \item{D}{
A matrix of distances obtained from an euclidean configuration
}
  \item{W}{
A matrix of weights
}
  \item{Model}{
Measurement level of the proximities
}
  \item{Standardize}{
Should the Disparities be standardized?
}
}
\details{
Calculation of disparities using standard  or monotone regression depending on the MDS model.
}
\value{
Returns the proximities.
}
\references{
Borg, I., & Groenen, P. J. (2005). Modern multidimensional scaling: Theory and applications. Springer.
}
\author{
Jose L. Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Function is used inside MDS or smacof
}