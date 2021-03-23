\name{plot.Ordinal.Logistic.Biplot}
\alias{plot.Ordinal.Logistic.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an ordinal Logistic Biplot
}
\description{
Plots an ordinal Logistic Biplot}
\usage{
\method{plot}{Ordinal.Logistic.Biplot}(x, A1 = 1, A2 = 2, 
ShowAxis = FALSE, margin = 0, PlotVars = TRUE, PlotInd = TRUE, 
LabelVars = TRUE, LabelInd = TRUE, mode = "a", CexInd = NULL, 
CexVar = NULL, ColorInd = NULL, ColorVar = NULL, SmartLabels = TRUE,
MinQualityVars = 0, dp = 0, PredPoints = 0, PlotAxis = FALSE, 
TypeScale = "Complete", ValuesScale = "Original", 
SizeQualInd = FALSE, SizeQualVars = FALSE, ColorQualInd = FALSE, 
ColorQualVars = FALSE, PchInd = NULL, PchVar = NULL, 
PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, 
ClustCenters = FALSE, UseClusterColors = TRUE, ClustLegend = TRUE,
ClustLegendPos = "topright", TextVarPos = 1, PlotSupVars = FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Plots and object of type "Ordinal.Logistic.Biplot"
}
  \item{A1}{
First dimension to plot
}
  \item{A2}{
Second dimension to plot
}
  \item{ShowAxis}{
Should the axis be shown
}
  \item{margin}{
Margin for the graph (in order to have space for the variable levels)
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{LabelVars}{
Should the variables be labelled?
}
  \item{LabelInd}{
Should the variables be labelled?
}
  \item{mode}{
Mode of the biplot (see the classical biplot)
}
  \item{CexInd}{
Type of marker used for the individuals
}
  \item{CexVar}{
Type of marker used for the variables
}
  \item{ColorInd}{
Colors used for the individuals
}
  \item{ColorVar}{
Colors used for the cariables
}
  \item{SmartLabels}{
Should smart placement for the labels be used?
}
  \item{MinQualityVars}{
Minimum quality of representation for a variable to be displayed
}
  \item{dp}{
Set of variables in which the individuals are projected
}
  \item{PredPoints}{
Set of points thet will be projected on all the variables
}
  \item{PlotAxis}{
Should the axis be plotted?
}
  \item{TypeScale}{
See continuous biplots
}
  \item{ValuesScale}{
See continuous biplots
}
  \item{SizeQualInd}{
Should the size of the labels and points be related to the quality of representation for individuals?
}
  \item{SizeQualVars}{
Should the size of the labels and points be related to the quality of representation for variables?
}
  \item{ColorQualInd}{
Should the intensity of the color of the labels and points be related to the quality of representation for individuals?
}
  \item{ColorQualVars}{
Should the intensity of the color of the labels and points be related to the quality of representation for variables?
}
  \item{PchInd}{
Markers for the individuals
}
  \item{PchVar}{
Markers for the individuals
}
  \item{PlotClus}{
Should the added clusters for the individuals be plotted?
}
  \item{TypeClus}{
Type of plot for the clusters. The types are "ch", "el" and "st" for "Convex Hull", "Ellipse" and "Star" repectively.
}
  \item{ClustConf}{
Confidence level for the cluster
}
  \item{ClustCenters}{
Should the centers of the clsters be plotted
}
  \item{UseClusterColors}{
Should the colors of the clusters be used to plot the individuals.
}
  \item{ClustLegend}{
Should a legend for the clusters be added?
}
  \item{ClustLegendPos}{
Position of the legend
}
  \item{TextVarPos}{
Position of the labels for the variables
}  
\item{PlotSupVars}{
Should the supplementary variables be plotted
}
  \item{\dots}{
Any other aditional parameters
}
}
\details{
Plots an ordinal Logistic Biplot
}
\value{
The plot ....
}
\references{
Vicente-Villardón, J. L., & Sánchez, J. C. H. (2014). Logistic Biplots for Ordinal Data with an Application to Job Satisfaction of Doctorate Degree Holders in Spain. arXiv preprint arXiv:1405.0294.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
    data(Doctors)
    olb = OrdLogBipEM(Doctors,dim = 2, nnodes = 10, initial=4,  tol = 0.001, 
    maxiter = 100, penalization = 0.1, show=TRUE)
    plot(olb, mode="s", ColorInd="gray", ColorVar=1:5)
}