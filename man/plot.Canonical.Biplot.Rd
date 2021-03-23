\name{plot.Canonical.Biplot}
\alias{plot.Canonical.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a Canonical Biplot
}
\description{
Plots a Canonical Biplot
}
\usage{
\method{plot}{Canonical.Biplot}(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups =
                    TRUE, PlotVars = TRUE, PlotInd = TRUE, WhatInds =
                    NULL, WhatVars = NULL, WhatGroups = NULL, IndLabels =
                    NULL, VarLabels = NULL, GroupLabels = NULL,
                    AbbreviateLabels = FALSE, LabelInd = TRUE, LabelVars =
                    TRUE, CexGroup = 1, PchGroup = 16, margin = 0.1,
                    AddLegend = FALSE, ShowAxes = FALSE, LabelAxes =
                    FALSE, LabelGroups = TRUE, PlotCircle = TRUE,
                    ConvexHulls = FALSE, TypeCircle = "M", ColorGroups =
                    NULL, ColorVars = NULL, LegendPos = "topright",
                    ColorInd = NULL, voronoi = TRUE, mode = "a", TypeScale
                    = "Complete", ValuesScale = "Original", MinQualityVars
                    = 0, dpg = 0, dpi = 0, dp = 0, PredPoints = 0,
                    PlotAxis = FALSE, CexInd = NULL, CexVar = NULL, PchInd
                    = NULL, PchVar = NULL, ColorVar = NULL, ShowAxis =
                    FALSE, VoronoiColor = "black", ShowBox = FALSE,
                    ShowTitle = TRUE, PlotClus = FALSE, TypeClus = "ch",
                    ClustConf = 1, ClustCenters = FALSE, UseClusterColors
                    = TRUE, CexClustCenters = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "Canonical.Biplot"
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
  \item{PlotGroups}{
Shoud the group centers be plotted?
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{WhatInds}{
Logical vector to control what individuals (Rows) are plotted. (Can be also a binary vector)
}
  \item{WhatVars}{
Logical vector to control what variables (Columns) are plotted. (Can be also a binary vector)
}
  \item{WhatGroups}{
Logical vector to control what groups are plotted. (Can be also a binary vector)
}
\item{IndLabels}{
A set of labels for the individuals. If NULL the default object labels are used
}
\item{VarLabels}{
A set of labels for the variables. If NULL the default object labels are used
}
\item{GroupLabels}{
A set of labels for the groups. If NULL the default object labels are used
}
  \item{AbbreviateLabels}{
Should labels be abbreviated?
}
  \item{LabelInd}{
Should the individuals be labeled?
}
  \item{LabelVars}{
Should the variables be labeled?
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
  \item{voronoi}{
Should the voronoi diagram with the prediction regións for each group be plotted?
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
A set of indices with the individuasl that will show the projections on the variables 
}
  \item{dp}{
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
  \item{VoronoiColor}{
Color for the Voronoi diagram
}
  \item{ShowBox}{
Should a box around the poitns be plotted?
} 
  \item{ShowTitle}{
Should the title be shown?
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
  \item{ClustCenters}{
Should the cluster centers be plotted?
}
 \item{UseClusterColors}{
Should the cluster colors be used in the plot
}
 \item{CexClustCenters}{
Size of the cluster centres
}
  \item{...}{
Any other graphical parameters 
}

}
\details{
The function plots the results of a Canononical Biplot. The coordinates for Groups, Individuals and Variables can be shown or not on the plot, each of the three can also be labeled separately.
The are parameters to control the way each different set of coordinates is plotted and labeled.

There are several modes for plotting the biplot.

"p".- Points (Rows and Columns are represented by points)

"a" .- Arrows (The traditional representation with points for rows and arrows for columns)

"b" .- The arrows for the columns are extended to both extremes of the plot and labeled outside the plot area.

"h" .- The arrows for the columns are extended to the positive extreme of the plot and labeled outside the plot area.

"ah" .- Same as arrows but labeled outside the plot area.

"s" .- The directions (or biplot axes) have a graded scale for prediction of the original values.

The \emph{TypeScale} argument applies only to the \emph{"s"} mode. There are three types:

"Complete" .- An equally spaced scale covering the whole range of the data is calculates.

"StdDev" .- Mean with one, two and three stadard deviations 

"BoxPlot" .- Box-Plot like Scale (Median, 25 and 75 percentiles, maximum and minimum values.)


The \emph{ValuesScale} argument applies only to the \emph{"s"} mode and controls if the labels show the  \emph{Original} ot \emph{Transformed} values.

Some of the initial transformations are not compatible with some of the types of biplots and scales. For example, It is not possible to recover by projection the original values when you double centre de data. In that case you have the residuals for interaction and only the transformed values make sense.
}
\value{
No value returned
}
\references{
Amaro, I. R., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Varas, M. J., Vicente-Tavera, S., Molina, E., & Vicente-Villardon, J. L. (2005). Role of canonical biplot method in the study of building stones: an example from Spanish monumental heritage. Environmetrics, 16(4), 405-419.

Santana, M. A., Romay, G., Matehus, J., Villardon, J. L., & Demey, J. R. (2009). simple and low-cost strategy for micropropagation of cassava (Manihot esculenta Crantz). African Journal of Biotechnology, 8(16).

}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
X=wine[,4:21]
canbip=CanonicalBiplot(X, group=wine$Group)
plot(canbip, TypeCircle="U")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Discriminant Analysis }
\keyword{ Canonical Variate Analysis}
\keyword{ MANOVA}
