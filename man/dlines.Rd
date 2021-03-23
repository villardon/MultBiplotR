\name{dlines}
\alias{dlines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Connects two sets of points by lines
}
\description{
Connects two sets of points by lines in a rowwise manner. Adapted from Graffelman(2013)
}
\usage{
dlines(SetA, SetB, lin = "dotted", color = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SetA}{
First set of points
}
  \item{SetB}{
Second set of points
}
  \item{lin}{
Line style.
}
  \item{color}{
Line color
}
  \item{...}{
Any other graphical parameters 
}
}
\details{
Connects two sets of points by lines 
}
\value{
NULL
}
\references{
Jan Graffelman (2013). calibrate: Calibration of Scatterplot and Biplot Axes. R package version 1.7.2.
  http://CRAN.R-project.org/package=calibrate
}
\author{
Based on Graffelman (2013)
}

\examples{
## No examples
}