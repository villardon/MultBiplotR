\name{plot.CanonicalDistanceAnalysis}
\alias{plot.CanonicalDistanceAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a Canonical Distance Analysis
}
\description{
Plots a Canonical Distance Analysis
}
\usage{
\method{plot}{CanonicalDistanceAnalysis}(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, 
ShowAxis = FALSE, ShowAxes = FALSE, LabelAxes = TRUE, margin = 0.1, 
PlotAxis = FALSE, ShowBox = TRUE, PlotGroups = TRUE, LabelGroups = TRUE, 
CexGroup = 1.5, PchGroup = 16, ColorGroup = NULL, voronoi = TRUE, 
VoronoiColor = "black", PlotInd = TRUE, LabelInd = TRUE, CexInd = 0.8, 
PchInd = 3, ColorInd = NULL, WhatInds = NULL, IndLabels = NULL, 
PlotVars = TRUE, LabelVar = TRUE, CexVar = NULL, PchVar = NULL, 
ColorVar = NULL, WhatVars = NULL, VarLabels = NULL, mode = "a", 
TypeScale = "Complete", ValuesScale = "Original", SmartLabels = TRUE, 
AddLegend = TRUE, LegendPos = "topright", PlotCircle = TRUE, 
ConvexHulls = FALSE, TypeCircle = "M", MinQualityVars = 0, dpg = 0, 
dpi = 0, PredPoints = 0, PlotClus = TRUE, TypeClus = "ch", ClustConf = 1, 
CexClustCenters = 1, ClustCenters = FALSE, UseClusterColors = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "CanonicalDistanceAnalysis"
}
  \item{A1}{
Dimension for the first axis. 1 is the default.
}
  \item{A2}{
Dimension for the second axis. 2 is the default.
}
  \item{ScaleGraph}{
Reescale the coordinates to optimal matching.
}
  \item{ShowAxis}{
Should the axis be shown?
}
  \item{ShowAxes}{
Not used
}
  \item{LabelAxes}{
Shoud the axis be labelled?
}
  \item{margin}{
Margin of the plot
}
  \item{PlotAxis}{
Should the axis be plotted?
}
  \item{ShowBox}{
Show a box around the plot
}
  \item{PlotGroups}{
Should the groups be plotted?
}
  \item{LabelGroups}{
Should the groups be labelled?
}
  \item{CexGroup}{
Sizes for the groups
}
  \item{PchGroup}{
Marks for the groups
}
  \item{ColorGroup}{
Colors for the groups
}
  \item{voronoi}{
Should a voronoi diagram separating the groups be plotted?
}
  \item{VoronoiColor}{
Color for the voronoi diagram
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelInd}{
Should the individuals be labelled?
}
  \item{CexInd}{
Sizes for the individuals
}
  \item{PchInd}{
Marks for the individuals
}
  \item{ColorInd}{
Colors for the individuals
}
  \item{WhatInds}{
What indivduals are plotted
}
  \item{IndLabels}{
Labels for the individuals
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{LabelVar}{
Should the variables be labelled?
}
  \item{CexVar}{
Sizes for the variables
}
  \item{PchVar}{
Marks for the variables
}
  \item{ColorVar}{
User colors for the variables. Default colors will be used if NULL.
}
  \item{WhatVars}{
What Variables are plotted
}
  \item{VarLabels}{
User labels for the variables
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
  \item{SmartLabels}{
%%     ~~Describe \code{SmartLabels} here~~
}
  \item{AddLegend}{
Should a legend be added?
}
  \item{LegendPos}{
Position of the legend
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
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be plotted
}
  \item{dpg}{
A set of indices with the variables that will show the projections of the gorups 
}
  \item{dpi}{
A set of indices with the individuasl that will show the projections on the variables 
}
\item{PredPoints}{
A vector with integers. The group centers listed in the vector are projected onto all the variables.
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
SIze of the cluster centers.
}
  \item{ClustCenters}{
Should the cluster centers be plotted?
}
 \item{UseClusterColors}{
Should the cluster colors be used in the plot
}
  \item{\dots}{
Any other graphical parameters 
}
}
\details{
Plots a Canonical Distance Analysis
}
\value{
The plot of a Canonical Distance Analysis
}
\references{
Gower, J. C. and Krzanowski, W. J. (1999). Analysis of distance for structured multivariate data and extensions to multivariate analysis of variance. Journal of the Royal Statistical Society: Series C (Applied Statistics), 48(4):505-519.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{plot.Canonical.Biplot}}
}
\examples{
# Not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
