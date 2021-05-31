\name{plot.PCoABootstrap}
\alias{plot.PCoABootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an object of class PCoABootstrap
}
\description{
Plots an object of class PCoABootstrap
}
\usage{
\method{plot}{PCoABootstrap}(x, F1=1, F2=2, Move2Center=TRUE, 
BootstrapPlot="Ellipse", confidence=0.95, Colors=NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "PCoABootstrap"
}
  \item{F1}{
First dimension to plot
}
  \item{F2}{
Second dimension to plot
}
  \item{Move2Center}{
Translate the ellipse center to the coordinates
}
  \item{BootstrapPlot}{
Type of Bootstrap plot to draw: "Ellipse", "ConvexHull", "Star"
}

  \item{confidence}{
Confidence level for the bootstrap plot
}
  \item{Colors}{
Colors of the objects
}
  \item{\dots}{
Additional parameters for graphical representations
}

}
\details{
Draws the bootstrap confidence regions for the coordinates of the points obtained from a Principal Coodinates Analysis
}
\value{
No value returned
}
\references{
J.R. Demey, J.L. Vicente-Villardon, M.P. Galindo, A.Y. Zambrano, Identifying molecular markers associated with classifications of genotypes by external logistic
biplot, Bioinformatics 24 (2008) 2832.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{data(spiders)
Dis=BinaryProximities(spiders)
pco=PrincipalCoordinates(Dis, Bootstrap=TRUE, BootstrapType="Products")
plot(pco, Bootstrap=TRUE)
}
}

