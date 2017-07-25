\name{plot.ContinuousBiplot}
\alias{plot.ContinuousBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a biplot for continuous data.
}
\description{
Plots a biplot for continuous data.
}
\usage{
\method{plot}{ContinuousBiplot}(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0,
                                  PlotVars = TRUE, PlotInd = TRUE, WhatInds = NULL,
                                  WhatVars = NULL, LabelVars = TRUE, LabelInd = TRUE,
                                  IndLabels = NULL, VarLabels = NULL, mode = "a", CexInd
                                  = NULL, CexVar = NULL, ColorInd = NULL, ColorVar =
                                    NULL, LabelPos = 1, SmartLabels = FALSE,
                                  MinQualityInds = 0, MinQualityVars = 0, dp = 0,
                                  PredPoints = 0, PlotAxis = FALSE, TypeScale =
                                    "Complete", ValuesScale = "Original", SizeQualInd =
                                    FALSE, SizeQualVars = FALSE, ColorQualInd = FALSE,
                                  ColorQualVars = FALSE, PchInd = NULL, PchVar = NULL,
                                  PlotClus = FALSE, TypeClus = "ch", ClustConf = 1,
                                  ClustLegend = FALSE, ClustLegendPos = "topright",
                                  ClustCenters = FALSE, UseClusterColors = TRUE, 
                                  PlotSupVars = FALSE, ShowBox=FALSE, nticks=5, NonSelectedGray=FALSE, 
                                  PlotUnitCircle=TRUE, PlotContribFA=TRUE, AddArrow=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
An object of class "Biplot"
}
  \item{A1}{
Dimension for the first axis. 1 is the default.
}
  \item{A2}{
Dimension for the second axis. 2 is the default.
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
}
  \item{LabelPos}{
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
  \item{ClustLegend}{
Should a legend for the clusters be plotted? Default FALSE
}
  \item{ClustLegendPos}{
Position of the legend for the clusters. Default "topright"
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
  \item{ShowBox}{
Should a box around the poitns be plotted?
} 
  \item{nticks}{
Number of ticks for the representation of the variables
} 
  \item{NonSelectedGray}{
The nonselected individuals and variables aplotted in light gray colors
} 
  \item{PlotUnitCircle}{
Plot the unit circle in the biplot for a Factor Analysis in which the lenght of the column arrows is smaller than 1 and is the quality of representation.
} 
  \item{PlotContribFA}{
Plot circles in the biplot for a Factor Analysis with different values of the quality of representation.
}
\item{AddArrow}{
Add an arrow to the representation of other modes of the biplot.
}
\item{...}{
Any other graphical parameters 
}
}
\details{
Plots a biplot for continuous data. The Biplot for continuous data is taken as the basis of the plot. If there are a mixture of different types of variables (binary, nominal, abundance, ...) are added to the biplot as supplementary parts.

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


It is possible to associate the color and the size of the points with the quality of representation. Bigger points correspond to better representation quality.

}
\value{
No value Returned
}
\references{

Gabriel, K. R. (1971). The biplot graphic display of matrices with application to principal component analysis. Biometrika, 58(3), 453-467.

Galindo Villardon, M. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Questiio. 1986, vol. 10, num. 1.

Vicente-Villardon, J. L., Galindo Villardon, M. P., & Blazquez Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.

Gower, J. C., & Hand, D. J. (1995). Biplots (Vol. 54). CRC Press.

Gower, J. C., Lubbe, S. G., & Le Roux, N. J. (2011). Understanding biplots. John Wiley & Sons.

Blasius, J., Eilers, P. H., & Gower, J. (2009). Better biplots. Computational Statistics & Data Analysis, 53(8), 3145-3158.
}
\author{
Jose Luis Vicente Villardon
}


\examples{
data(Protein)
bip=PCA.Biplot(Protein[,3:11])
## Biplot with scales on the variables
plot(bip, mode="s", margin=0.2, ShowAxis=FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
