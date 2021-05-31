\name{plot3dCanonicalBiplot}
\alias{plot3dCanonicalBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D Canonical Biplot
}
\description{
Plots a 3D Canonical Biplot
}
\usage{
plot3dCanonicalBiplot(Bip, A1 = 1, A2 = 2, A3 = 3, ScaleGraph = TRUE, 
PlotGroups = TRUE, PlotVars = TRUE, PlotInd = TRUE, LabelInd = TRUE, 
CexGroup = 1, PchGroup = 16, margin = 0.1, AddLegend = FALSE, 
ShowAxes = FALSE, LabelAxes = TRUE, LabelGroups = TRUE, PlotCircle = TRUE, 
ConvexHulls = FALSE, TypeCircle = "M", ColorGroups = NULL, ColorVars = NULL,
LegendPos = "topright", ColorInd = NULL, mode = "a", 
TypeScale = "Complete", ValuesScale = "Original", MinQualityVars = 0, 
dpg = 0, dpi = 0, PredPoints = 0, PlotAxis = FALSE, CexInd = NULL, 
CexVar = NULL, PchInd = NULL, PchVar = NULL, ColorVar = NULL, 
ShowAxis = TRUE, ColorAxis = "gray", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bip}{
An object of class "Canonical Biplot"
}
  \item{A1}{
Dimension for the first axis. 1 is the default.
}
  \item{A2}{
Dimension for the second axis. 2 is the default.
}
  \item{A3}{
Dimension for the third axis. 3 is the default.
}
  \item{ScaleGraph}{
Reescale the coordinates to optimal matching.
}
  \item{PlotGroups}{
Shoud the group centers be plotted?
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelInd}{
Should the individuals be labeled?
}
  \item{CexGroup}{
Sizes of the points for the groups
}
  \item{PchGroup}{
Markers for the group
}
  \item{margin}{
margin for the graph
}
  \item{AddLegend}{
Should a legend with the groups be added?
}
  \item{ShowAxes}{
Should outside axes be shown?
}
  \item{LabelAxes}{
Should outside axes be labelled?
}
  \item{LabelGroups}{
Should the groups be labeled?
}
  \item{PlotCircle}{
Should the confidence regions for the groups be plotted?
}
  \item{ConvexHulls}{
Should the convex hulls containing the individuals for each group be plotted?
}
  \item{TypeCircle}{
Type of confidence region: Univariate (U), Bonferroni(B), Multivariate (M) or Classical (C)
}
  \item{ColorGroups}{
User colors for the groups. Default colors will be used if NULL.
}
  \item{ColorVars}{
User colors for the variables. Default colors will be used if NULL.
}
  \item{LegendPos}{
Position of the legend.
}
  \item{ColorInd}{
User colors for the individuals. Default colors will be used if NULL.
}
 \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot"
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed"
}
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be plotted
}
  \item{dpg}{
A set of indices with the variables that will show the projections of the gorups 
}
  \item{dpi}{
A set of indices with the variables that will show the projections of the individuals
}
\item{PredPoints}{
A vector with integers. The group centers listed in the vector are projected onto all the variables.
}
  \item{PlotAxis}{
Not Used
}
  \item{CexInd}{
Size of the points for individuals.
}
  \item{CexVar}{
Size of the points for variables.
}
  \item{PchInd}{
Marhers of the points for individuals.
}
  \item{PchVar}{
Markers of the points for variables.
}
  \item{ColorVar}{
Colors of the points for variables.
}
  \item{ShowAxis}{
Should axis scales be shown?
}
  \item{ColorAxis}{
Color of the axis
}

  \item{\dots}{
Any aditional graphical parameters
}
}
\details{
The parameters are the same as in the 2D Canonical Biplot.
}
\value{
A 3D Canonical Biplot
}

\author{
Jose Luis Vicente Villardon
}
\seealso{
\code{\link{plot.Canonical.Biplot}}
}
\examples{
\donttest{data(wine)
X=wine[,4:21]
canbip=CanonicalBiplot(X, group=wine$Group)
plot3dCanonicalBiplot(canbip, TypeCircle="M")}
}