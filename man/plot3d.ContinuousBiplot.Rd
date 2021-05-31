\name{plot3d.ContinuousBiplot}
\alias{plot3d.ContinuousBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a classical biplot for continuous data
}
\description{
Plots a classical biplot for continuous data. 
}
\usage{
\method{plot3d}{ContinuousBiplot}(x, A1 = 1, A2 = 2, A3 = 3, 
ShowAxis = TRUE, margin = 0, PlotVars = TRUE, PlotInd = TRUE, 
WhatInds = NULL, WhatVars = NULL, LabelVars = TRUE, 
LabelInd = TRUE, IndLabels = NULL, VarLabels = NULL,
mode = "a", CexInd = NULL, CexVar = NULL, ColorInd = NULL, 
ColorVar = NULL, LabelPos = 1, SmartLabels = FALSE, 
MinQualityInds = 0, MinQualityVars = 0.3, dp = 0, 
PredPoints = 0, PlotAxis = FALSE, TypeScale = "Complete", 
ValuesScale = "Original", SizeQualInd = FALSE, 
SizeQualVars = FALSE, ColorQualInd = FALSE, 
ColorQualVars = FALSE, PchInd = NULL, PchVar = NULL, 
PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, 
ClustCenters = FALSE, UseClusterColors = TRUE, 
PlotSupVars = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "ContinuousBiplot""
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
  \item{ShowAxis}{
Logical variable to control if the coordinate axes should appear in the plot. The default value is FALSE because for most of the biplots its presence is irrelevant.
}
  \item{margin}{
Margin for the labels in some of the biplot modes (percentage of the plot width). Default is 0. Increase the value if the labels are not completely plotted.
}
  \item{PlotVars}{
Logical to control if the Variables (Columns) are plotted.
}
  \item{PlotInd}{
Logical to control if the Individuals (Rows) are plotted.
}
  \item{WhatInds}{
Logical vector to control what individuals (Rows) are plotted. (Can be also a binary vector)
}
  \item{WhatVars}{
Logical vector to control what variables (Columns) are plotted. (Can be also a binary vector)
}
  \item{LabelVars}{
Logical to control if the labels for the Variables are shown
}
  \item{LabelInd}{
Logical to control if the labels for the individuals are shown
}
\item{IndLabels}{
A set of labels for the individuals. If NULL the default object labels are used
}
\item{VarLabels}{
A set of labels for the variables. If NULL the default object labels are used
}
  \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{CexInd}{
Size for the symbols and labels of the individuals
}
  \item{CexVar}{
Size for the symbols and labels of the variables
}
  \item{ColorInd}{
Color for the symbols and labels of the individuals
}
  \item{ColorVar}{
Color for the symbols and labels of the variables
}  \item{LabelPos}{
Position of the labels in relation to the point. (Se the graphical parameter \code{pos} )
}
  \item{SmartLabels}{
Plot the labels in a smart way
}
  \item{MinQualityInds}{
Minimum quality of representation for an individual to be plotted
}
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be plotted
}
  \item{dp}{
A set of indices with the variables that will show the projections of the individuals
}
\item{PredPoints}{
A vector with integers. The row points listed in the vector are projected onto all the variables.
}
  \item{PlotAxis}{
Not Used
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot"
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed"
}
  \item{SizeQualInd}{
Should the size of the row points be related to their qualities of representation (predictiveness)?
}
  \item{SizeQualVars}{
Should the size of the column points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualInd}{
Should the color of the row points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualVars}{
Should the color of the column points be related to their qualities of representation (predictiveness)?
}
  \item{PchInd}{
Symbol for the row points. See \code{help(points)} for details.
}
  \item{PchVar}{
Symbol for the column points. See \code{help(points)} for details.
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
Should the cluster centers be plotted
}
 \item{UseClusterColors}{
Should the cluster colors be used in the plot
}
  \item{PlotSupVars}{
Should the supplementary variables be plotted?
} 
  \item{\dots}{
Any other graphical parameters 
}
}
\details{
The parameters are the same as the ones  for the 2D biblot.
}
\value{
A 3D Biplot
}

\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
data(Protein)
bip=PCA.Biplot(Protein[,3:11])
plot3d.ContinuousBiplot(bip, mode="s", margin=0.2, ShowAxis=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ 3D Biplot }
