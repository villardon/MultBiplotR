\name{plot.Binary.Logistic.Biplot}
\alias{plot.Binary.Logistic.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the results of a Binary Logistic  Biplot
}
\description{
Plots the results of a Binary Logistic  Biplot
}
\usage{
\method{plot}{Binary.Logistic.Biplot}(x, F1 = 1, F2 = 2, ShowAxis = FALSE, margin = 0, 
PlotVars = TRUE, PlotInd = TRUE, WhatRows = NULL, WhatCols = NULL, 
LabelRows = TRUE, LabelCols = TRUE, ShowBox = FALSE, RowLabels = NULL, 
ColLabels = NULL, RowColors = NULL, ColColors = NULL, Mode = "s", 
TickLength = 0.01, RowCex = 0.8, ColCex = 0.8, SmartLabels = FALSE, 
MinQualityRows = 0, MinQualityCols = 0, dp = 0, PredPoints = 0, 
SizeQualRows = FALSE, SizeQualCols = FALSE, ColorQualRows = FALSE, 
ColorQualCols = FALSE, PchRows = NULL, PchCols = NULL, PlotClus = FALSE, 
TypeClus = "ch", ClustConf = 1, Significant = TRUE, alpha = 0.05, 
Bonferroni = TRUE, PlotSupVars = TRUE, AbbreviateLabels = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class Binary.Logistic.Biplot
}
  \item{F1}{
Dimension for the first axis of the representation. Default = 1
}
  \item{F2}{
Dimension for the second axis of the representation. Default = 2
}
  \item{ShowAxis}{
Should the axis of the representation be shown?
}
  \item{margin}{
Margin of the plot as a percentage. It gets some space for the labels.
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{WhatRows}{
What Rows should be plotted. A binary vector containing which rows (individuals) should be plotted (1) and which should not (0).
}
  \item{WhatCols}{
What Columns should be plotted. A binary vector containing which columns (variables) should be plotted (1) and which should not (0).
}
  \item{LabelRows}{
Should the individuals be labeled?
}
  \item{LabelCols}{
Should the individuals be labeled?
}
  \item{ShowBox}{
Should a box around the points be plotted?
}
  \item{RowLabels}{
A vector of row labels. If NULL the labels contained in the object will be used.
}
  \item{ColLabels}{
A vector of column labels. If NULL the labels contained in the object will be used.
}
  \item{RowColors}{
A vector of alternative row colors. 
}
  \item{ColColors}{
A vector of alternative column colors. 
}
  \item{Mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{TickLength}{
Length of the scale ticks for the biplot variables.
}
  \item{RowCex}{
Cex (Size) of the rows (marks and labels). Can be a single common size for all the points or a vector with individual sizes.
}
  \item{ColCex}{
Cex (Size) of the columns (marks and labels). Can be a single common size for all the points or a vector with individual sizes.
}
  \item{SmartLabels}{
Should the labels be placed in a smart way?
}
  \item{MinQualityRows}{
Minimum quality of the rows to be plotted. (Between 0 and 1)
}
  \item{MinQualityCols}{
Minimum quality of the columns to be plotted. (Between 0 and 1)
}
  \item{dp}{
A vector of variable indices to project all the individuals onto each variable of the vector.
}
  \item{PredPoints}{
A vector of row indices to project onto each variable.
}
  \item{SizeQualRows}{
Should the size of the Row points be related to its quality?
}
  \item{SizeQualCols}{
Should the size of the Column points be related to its quality?
}
  \item{ColorQualRows}{
Should the color of the Row points be related to its quality?
}
  \item{ColorQualCols}{
Should the color of the Column points be related to its quality?
}
  \item{PchRows}{
Marks for the rows (numbers). Can be a single common mark for all the points or a vector with individual marks.
}
  \item{PchCols}{
Marks for the columns (numbers). Can be a single common mark for all the points or a vector with individual marks.
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
  \item{Significant}{
Should only the significant variables be plotted?
}
  \item{alpha}{
Signification level.
}
  \item{Bonferroni}{
Should the Bonferroni correction be used?
}
  \item{PlotSupVars}{
Should the Supplementary variables be plotted?
}
  \item{AbbreviateLabels}{
Should labels be abbreviated?
}
  \item{\dots}{
Any other graphical parameter.
}
}
\details{
Plots a biplot for binary data. The Biplot for binary data is taken as the basis of the plot. If there are a mixture of different types of variables (binary, nominal, abundance, ...) are added to the biplot as supplementary parts.

There are several modes for plotting the biplot.
"p".- Points (Rows and Columns are represented by points)

"a" .- Arrows (The traditional representation with points for rows and arrows for columns)

"b" .- The arrows for the columns are extended to both extremes of the plot and labeled outside the plot area.

"h" .- The arrows for the columns are extended to the positive extreme of the plot and labeled outside the plot area.

"ah" .- Same as arrows but labeled outside the plot area.

"s" .- The directions (or biplot axes) have a graded scale for prediction of the original values.
}
\value{
The plot of the biplot.
}
\references{
Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Análisis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(spiders)
X=Dataframe2BinaryMatrix(spiders)

logbip=BinaryLogBiplotGD(X,penalization=0.1)
plot(logbip, Mode="a")
summary(logbip)
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
