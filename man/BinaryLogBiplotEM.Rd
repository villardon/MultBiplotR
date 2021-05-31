\name{BinaryLogBiplotEM}
\alias{BinaryLogBiplotEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary logistic biplot with the EM algorithm.
}
\description{
Binary logistic biplot with the EM algorithm
}
\usage{
BinaryLogBiplotEM(x, freq = matrix(1, nrow(x), 1), aini = NULL,
dimens = 2, nnodos = 15, tol = 1e-04, maxiter = 100, penalization = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A binary data matrix
}
  \item{freq}{
A vector of frequencies.
}
  \item{aini}{
Initial values for the row coordinates.
}
  \item{dimens}{
Dimension of the solution.
}
  \item{nnodos}{
Number of nodes for the gaussian quadrature
}
  \item{tol}{
Tolerance
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{penalization}{
Penalization for the fit (ridge)
}
}
\details{
Binary logistic biplot with the EM algorithm based on marginal maximum likelihood.
}
\value{
A logistic biplot object.
}
\references{
Vicente-Villardón, J. L., Galindo-Villardón, M. P., & Blázquez-Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.
}
\author{
Jose Luis Vicente-Villardon
}

\examples{
# Not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
