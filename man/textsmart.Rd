\name{textsmart}
\alias{textsmart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Labels of a Scatter}
\description{
Plots labels of points in a scattergram. labels for points with positive \code{x} are placed on the right of 
the points, and labels for points with negative values on the left. 
}
\usage{
textsmart(A, Labels, CexPoints, ColorPoints, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
Coordinates of the points for the scaterrgram
}
  \item{Labels}{
Labels for the points
}
  \item{CexPoints}{
Size of the labels
}
  \item{ColorPoints}{
Colors of the labels
}
  \item{...}{
Aditional graphical arguments
}
}
\details{
The function is used to improve the readability of the labels in a scatergram.
}
\value{
No value returned
}

\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{plot.Principal.Coordinates}}
}
\examples{
data(spiders)
dist=BinaryProximities(spiders)
pco=PrincipalCoordinates(dist)
plot(pco, SmartLabels =TRUE)
}

