\name{BinaryLogBiplotGD}
\alias{BinaryLogBiplotGD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Binary Logistic Biplot with Gradient Descent Estimation
}
\description{
Binary Logistic Biplot with Gradient Descent Estimation. An external optimization function is used to calculate the parameters.
}
\usage{
BinaryLogBiplotGD(X, freq = matrix(1, nrow(X), 1), dim = 2, tolerance = 1e-04, penalization = 0.01, RotVarimax = FALSE, seed = NULL, OptimMethod = "CG", Initial = "random", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A binary data matrix
}
  \item{freq}{
Frequencies of each row. When adequate.}
  \item{dim}{
Dimension of the final solution.
}
  \item{tolerance}{
Tolerance for convergence of the algorithm.
}
  \item{penalization}{
Ridge penalization constant.
}
  \item{RotVarimax}{
Should the final solution be rotated.
}
  \item{seed}{
Seed for the random numbers. Used for reproductibility.
}
  \item{OptimMethod}{
Optimization method used by \code{optimr}.
}
  \item{Initial}{
Initial configuration to start the iterations.
}
  \item{\dots}{
Aditional parameters used by the optimization function.
}
}
\details{
Fits a binary logistic biplot using gradient descent. The general function \code{optimr} is used to optimize the loss function. Conjugate gradien is used as a default although other alternatives can be USED.
}
\value{
An object of class "Binary.Logistic.Biplot".
}
\references{
Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Análisis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.
}
\author{
Jose Luis Vicente-Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(spiders)
X=Dataframe2BinaryMatrix(spiders)

logbip=BinaryLogBiplotGD(X,penalization=0.1)
plot(logbip, Mode="a")
summary(logbip)
}
