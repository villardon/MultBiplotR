\name{AddBinVars2Biplot}
\alias{AddBinVars2Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add suplementary binary variables to a biplot
}
\description{
Add suplementary binary variables to a biplot of any kind
}
\usage{
AddBinVars2Biplot(bip, Y, IncludeConst = TRUE, penalization = 0.2, 
freq = NULL, tolerance = 1e-05, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
A biplot object
}
  \item{Y}{
Matrix of binary variables to add
}
  \item{IncludeConst}{
Should include a constant in the fit
}
  \item{penalization}{
Penalization for the fit
}
  \item{freq}{
frequencies for each row of Y. By default is 1.
}
  \item{tolerance}{
Tolerance for the fit
}
  \item{maxiter}{
Maximum number of iterations
}
}

\details{
Fits binary variables to an existing biplot using penalized logistic regression.
}
\value{
The biplot object with supplementary binary variables added.
}
\references{
Vicente-Villardón, J. L., & Hernández-Sánchez, J. C. (2020). External Logistic Biplots for Mixed Types of Data. In Advanced Studies in Classification and Data Science (pp. 169-183). Springer, Singapore.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
## No examples yet
}
