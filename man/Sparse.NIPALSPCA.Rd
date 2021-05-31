\name{Sparse.NIPALSPCA}
\alias{Sparse.NIPALSPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse version of the NIPALS algorithm for PCA.
}
\description{
Sparse version of the NIPALS algorithm for PCA.
}
\usage{
Sparse.NIPALSPCA(X, dimens = 2, tol = 1e-06, maxiter = 1000, lambda = 0.02)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The data matrix.
}
  \item{dimens}{
The dimension of the solution
}
  \item{tol}{
Tolerance of the algorithm.
}
  \item{maxiter}{
Maximum number of iteratios.
}
  \item{lambda}{
Value used for sparsity
}
}
\details{
Sparse version of the NIPALS algorithm for the singular value decomposition that allows for the construction of PCA and Biplot.
}
\value{
The singular value decomposition

\item{u}{The coordinates of the rows (standardized)}
\item{d }{The singuklar values}
\item{v }{The coordinates of the columns (standardized)}
}
\references{
Have to be written
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
