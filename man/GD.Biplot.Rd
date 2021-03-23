\name{GD.Biplot}
\alias{GD.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot for continuous data based on gradient descent methods
}
\description{
Biplot for continuous data based on gradient descent methods.
}
\usage{
GD.Biplot(X, dimension = 2, Scaling = 5, 
         lambda = 0.01, OptimMethod = "CG", 
         Orthogonalize = FALSE, Algorithm = "Alternated", 
         sup.rows = NULL, sup.cols = NULL,
         grouping = NULL, tolerance = 1e-04, 
         num_max_iters = 300, Initial = "random")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix with continuous variables.
}
  \item{dimension}{
Dimension of the final solution.
}
  \item{Scaling}{
Transformation of the raw data matrix before the calculation of the biplot.
}
  \item{lambda}{
  Constant for the ridge Penalization
}
  \item{OptimMethod}{
Optimization method passed to the \code{optimr} function. By default is CG (Conjugate Gradient).
}
  \item{Orthogonalize}{
Should the solution be ortogonalized.
}
  \item{Algorithm}{
Algorithm to calculate the Biplot. (Alternated, Joint, Recursive)
}
  \item{sup.rows}{
Supplementary Rows. (not working now)
}
  \item{sup.cols}{
Supplementary Columns. (not working now)
}
  \item{grouping}{
Grouping factor for the within groups transformation.
}
  \item{tolerance}{
Tolerance for convergence
}
  \item{num_max_iters}{
Maximum number of iterations.
}
  \item{Initial}{
Initial Configuration 
}
}
\details{
The function calculates a bilot using gradient descent methods. The function \code{optimr} is used to optimize the loss function. By default CG (Conjugate Gradient) method is used althoug other possibilities can be used.
}
\value{
An object of class "ContinuousBiplot" is returned.
}

\author{
Jose Luis Vicente Villardon
}


\examples{
data("Protein")
X=Protein[,3:11]
gdbip=GD.Biplot(X, dimension=2, Algorithm="Joint", 
Orthogonalize=FALSE, lambda=0.3, Initial="random")
plot(gdbip)
summary(gdbip)
}

