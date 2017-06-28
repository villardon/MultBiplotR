\name{plot.MGC}
\alias{plot.MGC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the results of Model-Based Gaussian Clustering algorithms
}
\description{
PLots an object of type MGC (Model-based Gaussian Clustering)
}
\usage{
\method{plot}{MGC}(x, vars = NULL, groups = x$Classification, CexPoints = 0.2, Confidence = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of type MGC
}
  \item{vars}{
A subset of indices of the variables to be plotted
}
  \item{groups}{
A factor containing groups to represent. Usually the clusters obtained from the algorithm.
}
  \item{CexPoints}{
Size of the points.
}
  \item{Confidence}{
Confidence of the ellipses
}
  \item{\dots}{
Anay additional graphical parameters
}
}
\details{
PLots an object of type MGC (Model-based Gaussian Clustering) using a splom plot.
}
\value{
No value returned
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(iris)
}
