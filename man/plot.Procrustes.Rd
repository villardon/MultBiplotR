\name{plot.Procrustes}
\alias{plot.Procrustes}
\title{Plots an object of class "Procrustes"}

\description{
  Plots Simple Procrustes Analysis
}

\usage{
\method{plot}{Procrustes}(x, F1=1, F2=2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "Procrustes"
}
  \item{F1}{
First dimenssion of the plot
}
  \item{F2}{
Second dimenssion of the plot
}

  \item{\dots}{
Additional parameters for graphical representations
}
}
\details{
Graphical representation of an Orthogonal Procrustes Analysis. 
}
\value{
No value is returned
}

\author{
Jose Luis Vicente-Villardon
}


\seealso{
\code{\link{BinaryProximities}}
}
\examples{
data(spiders)
dist=BinaryProximities(spiders)
pco=PrincipalCoordinates(dist)
plot(pco)
}
