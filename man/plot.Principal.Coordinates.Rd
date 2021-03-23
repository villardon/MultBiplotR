\name{plot.Principal.Coordinates}
\alias{plot.Principal.Coordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plots an object of class Principal.Coordinates}
\description{
Plots an object of class Principal.Coordinates
}
\usage{
\method{plot}{Principal.Coordinates}(x, A1 = 1, A2 = 2, LabelRows = TRUE, 
WhatRows = NULL, RowCex = 1, RowPch = 16, Title = "", RowLabels = NULL, 
RowColors = NULL, ColColors = NULL, ColLabels = NULL, SizeQualInd = FALSE, 
SmartLabels = TRUE, ColorQualInd = FALSE, ColorQual = "black", PlotSup = TRUE, 
Bootstrap = FALSE, BootstrapPlot = c("Ellipse", "CovexHull", "Star"), 
margin = 0, PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, 
CexClustCenters = 1, LegendClust = TRUE, ClustCenters = FALSE, 
UseClusterColors = TRUE, ShowAxis = FALSE, PlotBinaryMeans = FALSE, 
MinIncidence = 0, ShowBox = FALSE, ColorSupContVars = NULL, 
ColorSupBinVars = NULL, ColorSupOrdVars = NULL, TypeScale = "Complete", 
SupMode = "s", PlotSupVars = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class "Principal.Coordinates"
}
  \item{A1}{
First dimenssion of the plot
}
  \item{A2}{
Second dimenssion of the plot
}
  \item{LabelRows}{
Controls if the points are labelled. Usually TRUE.
}
  \item{WhatRows}{
What Rows to plot. A vector of 0/1 elements. If NULL all rows are plotted
}
  \item{RowCex}{
Size of the points. Can be a single number or a vector.
}
  \item{RowPch}{
Symbols for the points.
}
  \item{Title}{
Title for the graph
}
  \item{RowLabels}{
Labels for the rows. If NULL row names of the data matrix are used.
}
  \item{RowColors}{
Colors for the rows. If NULL row deafault colors are assigned. Can be a single value or avector of colors.
}
  \item{ColColors}{
Colors for the columns (Variables)
}
  \item{ColLabels}{
Labels for the columns (Variables)
}
  \item{SizeQualInd}{
Controls if the size of points depends on the quality of representation.
}
  \item{SmartLabels}{
Controls the way labels are plotted on the graph. If TRUE labels for points with positive x values
are placed to the right of the point and labels for points with negative values to the left
}
  \item{ColorQualInd}{
Controls if the color of the points depends on the quality of representation.
}
  \item{ColorQual}{
Darker color for the quality scale.
}
  \item{PlotSup}{
Controls if the supplementary points are plotted.
}
  \item{Bootstrap}{
Controls if the bootstrap points are plotted.
}
  \item{BootstrapPlot}{
Type of plot of the Bootstrap Information. The types are "Ellipse", "CovexHull" or "Star".
}
  \item{margin}{
Margin for the graph.
}
  \item{PlotClus}{
Should the clusters be plotted?
}
  \item{TypeClus}{
Type of plot for the clusters. ("ch"- Convex Hull, "el"- Ellipse or "st"- Star)
}
  \item{ClustConf}{
Percent of points included in the cluster. only the ClusConf percent of the points nearest to the center will be used to calculate the cluster
}
  \item{CexClustCenters}{
Size of the cluster centers
}
  \item{LegendClust}{
Legends for the clusters
}
 \item{ClustCenters}{
Should the cluster centers be plotted
}
 \item{UseClusterColors}{
Should the cluster colors be used in the plot
}
  \item{ShowAxis}{
%%     ~~Describe \code{ShowAxis} here~~
}
  \item{PlotBinaryMeans}{
%%     ~~Describe \code{PlotBinaryMeans} here~~
}
  \item{MinIncidence}{
%%     ~~Describe \code{MinIncidence} here~~
}
  \item{ShowBox}{
%%     ~~Describe \code{ShowBox} here~~
}
  \item{ColorSupContVars}{
%%     ~~Describe \code{ColorSupContVars} here~~
}
  \item{ColorSupBinVars}{
%%     ~~Describe \code{ColorSupBinVars} here~~
}
  \item{ColorSupOrdVars}{
%%     ~~Describe \code{ColorSupOrdVars} here~~
}
  \item{TypeScale}{
%%     ~~Describe \code{TypeScale} here~~
}
  \item{SupMode}{
%%     ~~Describe \code{SupMode} here~~
}
  \item{PlotSupVars}{
%%     ~~Describe \code{PlotSupVars} here~~
}
  \item{\dots}{
Additional parameters for graphical representations
}
}
\details{
Graphical representation of an Principal coordinates Analysis controlling visual aspects of the plot as colors, symbols or sizes of the points. 
}
\value{
No value is returned
}
\references{
J.R. Demey, J.L. Vicente-Villardon, M.P. Galindo, A.Y. Zambrano, Identifying molecular markers associated with classifications of genotypes by external logistic
biplot, Bioinformatics 24 (2008) 2832.
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

