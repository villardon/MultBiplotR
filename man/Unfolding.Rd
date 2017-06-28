\name{Unfolding}
\alias{Unfolding}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Unfolding para vegetacion
}
\description{
Unfolding para vegetacion
}
\usage{
Unfolding(A, ENV = NULL, TransAbund = "Gaussian", offset = 0.5, weight = "All_1", Constrained = FALSE, TransEnv = "Standardize columns", InitConfig = "SVD", model = "Ratio", condition = "Columns", r = 2, maxiter = 100, tolerance = 1e-05, lambda = 1, omega = 0, plot = FALSE)
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
%%     ~~Describe \code{weight} here~~
}
  \item{Constrained}{
%%     ~~Describe \code{Constrained} here~~
}
  \item{TransEnv}{
%%     ~~Describe \code{TransEnv} here~~
}
  \item{InitConfig}{
%%     ~~Describe \code{InitConfig} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{condition}{
%%     ~~Describe \code{condition} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{tolerance}{
%%     ~~Describe \code{tolerance} here~~
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
%%  ~~ If necessary, more details than the description above ~~
}
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
data("SpidersSp")
unf=Unfolding(SpidersSp)
plot(unf)
Genefold(SpidersSp)
}

