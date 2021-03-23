\name{DensityBiplot}
\alias{DensityBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds Non-parametric densities to a biplot. Separated densities are calculated for different clusters
}
\description{
Adds Non-parametric densities to a biplot. Separated densities are calculated for different clusters
}
\usage{
DensityBiplot(X, y = NULL, grouplabels = NULL, ncontours = 6, 
groupcolors = NULL, ncolors=20, ColorType=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Two dimensional coordinates of the points in a biplot (or any other)
}
  \item{y}{
A factor containing clusters or groups for separate densities. 
}
  \item{grouplabels}{
Labels for the groups
}
  \item{ncontours}{
Number of contours to represent on the biplot
}
  \item{groupcolors}{
Colors for the groups
}
  \item{ncolors}{
Number of colors for the density patterns
}
  \item{ColorType}{
One of the following: "1" = rainbow, "2" = heat.colors, "3" = terrain.colors, "4" = topo.colors, "5" = cm.colors
}      
}
\details{
Non parametric densities are used to investigate the concentration of row points on different areas of the biplot representation. The densities can be calculated for different groups or clusters in order to investigate if individuals with differnt characteristics are concentrated on particular areas of the biplot. The procedure is particularly useful with a high number of individuals.
}
\value{
No value returned. It has effect on the graph.
}
\references{

Gower, J. C., Lubbe, S. G., & Le Roux, N. J. (2011). Understanding biplots. John Wiley & Sons.
}
\author{
Jose Luis Vicente Villardon
}


\examples{
bip=PCA.Biplot(iris[,1:4])
plot(bip, mode="s", CexInd=0.1)
% DensityBiplot(X=bip$RowCoordinates[,1:2], y=iris[,5], groupcolors=c("green", "blue", "red"))

}

