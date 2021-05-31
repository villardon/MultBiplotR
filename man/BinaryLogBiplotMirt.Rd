\name{BinaryLogBiplotMirt}
\alias{BinaryLogBiplotMirt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary logistic biplot with Item Response Theory.
}
\description{
Binary logistic biplot with Item Response Theory.
}
\usage{
BinaryLogBiplotMirt(x, dimens = 2, tolerance = 1e-04, 
maxiter = 30, penalization = 0.2, Rotation = "varimax", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The binary Data matrix
}
  \item{dimens}{
Dimension of the solution
}
  \item{tolerance}{
Tolerance of the algorithm
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{penalization}{
Rige Penalization
}
  \item{Rotation}{
Should a rotation be applied?
}
  \item{\dots}{
Aditional argumaents.
}
}
\details{
Binary logistic biplot with Item Response Theory.
}
\value{
A logistic biplot object.
}
\references{
Vicente-Villardón, J. L., Galindo-Villardón, M. P., & Blázquez-Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.
}
\author{
Jose Luis Vicente Villardon
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
