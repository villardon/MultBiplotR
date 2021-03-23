\name{Unfolding}
\alias{Unfolding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional Unfolding
}
\description{
Multidimensional Unfolding with some adaptations for vegetation analysis
}
\usage{
Unfolding(A, ENV = NULL, TransAbund = "Gaussian Columns", offset = 0.5, 
weight = "All_1", Constrained = FALSE, 
TransEnv = "Standardize columns", 
InitConfig = "SVD", model = "Ratio", 
condition = "Columns", Algorithm = "SMACOF", 
OptimMethod = "CG", r = 2, maxiter = 100, 
tolerance = 1e-05, lambda = 1, omega = 0, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
The original proximities matrix
}
  \item{ENV}{
The matrix of environmental variables
}
  \item{TransAbund}{
Initial transformation of the abundances : "None", "Gaussian", "Column Percent", "Gaussian Columns", "Inverse Square Root", "Divide by Column Maximum")
}
  \item{offset}{
offset is the quantity added to the zeros of the table
}
  \item{weight}{
A matrix of weights for each cell of the table
}
  \item{Constrained}{
Should fit a constrained analysis
}
  \item{TransEnv}{
Transformation of the environmental variables
}
  \item{InitConfig}{
Init configuration for the algorithm
}
  \item{model}{
Type of model to be fitted: "Identity", "Ratio", "Interval" or "Ordinal".
}
  \item{condition}{
"Matrix", "Columns" to condition to the whole matrix or to each column
}
  \item{Algorithm}{
Algorithm to fit the model: "SMACOF", "GD", "Genefold"
}
  \item{OptimMethod}{
Optimization method for gradient descent
}
  \item{r}{
Dimension of the solution
}
  \item{maxiter}{
Maximum number of iterations in the algorithm
}
  \item{tolerance}{
Tolerace for the algorithm
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{omega}{
%%     ~~Describe \code{omega} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
}
\details{
ological data}
\value{
An object of class "Unfolding"
}
\references{
Ver Articulos
}
\author{
Jose Luis Vicente Villardon
}

\examples{
unf=Unfolding(SpidersSp, ENV=SpidersEnv, model="Ratio", Constrained = FALSE, condition="Matrix")
plot(unf, PlotTol=TRUE, PlotEnv = FALSE)
plot(unf, PlotTol=TRUE, PlotEnv = TRUE)
cbind(unf$QualityVars, unf$Var_Fit)
unf2=Unfolding(SpidersSp, ENV=SpidersEnv, model="Ratio", Constrained = TRUE, condition="Matrix")
plot(unf2, PlotTol=FALSE, PlotEnv = TRUE, mode="s")
cbind(unf2$QualityVars, unf2$Var_Fit)
}

