\name{plot.External.Binary.Logistic.Biplot}
\alias{plot.External.Binary.Logistic.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an External Logistic Biplot for binary data
}
\description{
Plot of an External Binary Logistic Biplot with many arguments controling different aspects of the representation
}
\usage{
\method{plot}{External.Binary.Logistic.Biplot}(x, F1 = 1, F2 = 2, 
                    ShowAxis = FALSE, margin = 0.1,
                    PlotVars = TRUE, PlotInd = TRUE, WhatRows = NULL,
                    WhatCols = NULL, LabelRows = TRUE, LabelCols = TRUE,
                    RowLabels = NULL, ColLabels = NULL, RowColors = NULL,
                    ColColors = NULL, Mode = "s", TickLength = 0.01,
                    RowCex = 0.8, ColCex = 0.8, SmartLabels = FALSE,
                    MinQualityRows = 0, MinQualityCols = 0, dp = 0,
                    PredPoints = 0, SizeQualRows = FALSE, ShowBox = FALSE,
                    SizeQualCols = FALSE, ColorQualRows = FALSE,
                    ColorQualCols = FALSE, PchRows = NULL, PchCols = NULL,
                    PlotClus = FALSE, TypeClus = "ch", ClustConf = 1,
                    Significant = FALSE, alpha = 0.05, Bonferroni = FALSE,
                    PlotSupVars = TRUE, ...)
                    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of type \code{External.Binary.Logistic.Biplot} 
}
  \item{F1}{
Latent factor to represent at the X axis 
}
  \item{F2}{
Latent factor to represent at the Y axis 
}
  \item{ShowAxis}{
Should the axis be plotted? 
}
  \item{margin}{
Margin for the labels in some of the biplot modes (percentage of the plot width). Default is 0. Increase the value if the labels are not completely plotted.
}
  \item{PlotVars}{
Should Variables be plotted
}
  \item{PlotInd}{
Should Individuals be plotted 
}
  \item{WhatRows}{
A binary vector (0 and 1) that indicates if each individual row should be plotted or not 
}
  \item{WhatCols}{
A binary vector (0 and 1) that indicates if each individual column should be plotted or not 
}
  \item{LabelRows}{
Should Variables be labelled
}
  \item{LabelCols}{
Should Individuals be labelled 
}
  \item{RowLabels}{
A vector of Labels for the rows if you do not want to use the data labels
}
\item{ColLabels}{
A vector of Labels for the columns if you do not want to use the data labels
}
\item{RowColors}{
A vector of colors for the rows 
}
\item{ColColors}{
A vector of colors for the rows 
}
  \item{Mode}{
Mode of the biplot: "p", "a", "b", "ah" and "s". See details.
}
  \item{TickLength}{
Lenght of the tick marks. Depends on the scale of the graph.
}
\item{RowCex}{
A scalar or a vector containing the sizes of the poitns ans labels for the rows. Default value is 0.8 if the sizes are not provided.
}
\item{ColCex}{
A scalar or a vector containing the sizes of the poitns ans labels for the columns. Default value is 0.8 if the sizes are not provided.
}
\item{SmartLabels}{
Plot the labels in a smart way
}
\item{MinQualityRows}{
Minimum quality of representation for a row or individual to be plotted
}
\item{MinQualityCols}{
Minimum quality of representation for a column or variable to be plotted
}
\item{dp}{
"Drop Points" on the variables, a vector with integers. The row points are projected on the directions of the variables listed in the vector.
}
\item{PredPoints}{
A vector with integers. The row points listed in the vector are projected onto all the variables.
}

  \item{SizeQualRows}{
Should the size of the row points be related to their qualities of representation (predictiveness)?
}
  \item{ShowBox}{
Should abox around the point be displayed?
}
  \item{SizeQualCols}{
Should the size of the column points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualRows}{
Should the color of the row points be related to their qualities of representation (predictiveness)?
}
  \item{ColorQualCols}{
Should the color of the column points be related to their qualities of representation (predictiveness)?
}
 \item{PchRows}{
Symbol for the row points. See \code{help(points)} for details.
}
  \item{PchCols}{
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
  \item{Significant}{
If TRUE, only the significant variables are plotted 
}
  \item{alpha}{
Significance Level
}
  \item{Bonferroni}{
Should the Bonferroni correction be used
}   \item{PlotSupVars}{
Should supplementary variables be plotted
}    
\item{...}{
Any other graphical parameter you want to use
}
}
\details{
The logistic regression equation predicts the probability that a caracter will be present in an individual. Geometrically the y´s can be represented as point in the reduced dimension space and the b's are the vectors showing the directions that best predict the probability of presence of each allele . For a com-plete explanation of the geometrical properties of the ELB see Vicente-Villardón et al (2006).
The prediction of the probabilities is made in the same way as in a linear Biplot, i. e., the projection of a genotype point on the direction of an variable vector predicts the probability of presence of that variable in the individual. To facilitate the interpretation of the graph, fixed prediction probabilities points are situated on each allele vector. To simplify the graph, in our ap-plication, a vector joining the points for 0.5 and 0.75 are placed; this shows the cut point for prediction of presence and the direction of increasing probabilities. The length of the vector can be interpreted as an inverse measure of the discriminatory power of the alleles or bands, in the sense that shorter vectors correspond to alleles that better differentiate individuals. Two alleles pointing in the same direction are highly correlated, two alleles pointing in opposite directions are negatively correlated, and two alleles forming an angle close to 90º are not correlated. A more complete scale with probabilities from 0.1 to 0.9 can also be plotted with this function.
For each variable, the ordination diagram can be divided into two separate regions predicting presence or absence, the two regions are separated by the line that is perpendicular to the variable vector in the Biplot and cuts the vector in the point predicting 0.5.  The variables associated to the configuration are those that predict the presences adequately.
In a practical situation not all the variables are associated to the ordination. Due to the high number usually studied, it is convenient to situate on the graph only those that are related to the configuration, i. e., those that have an adequate goodness of fit after adjusting the logistic regression.

}
\value{
No value returned
}
\references{
Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.

Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Analysis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.
}
\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ExternalBinaryLogisticBiplot}}
}
\examples{
data(spiders)
dist=BinaryProximities(spiders)
pco=PrincipalCoordinates(dist)
pcobip=ExternalBinaryLogisticBiplot(pco)
plot(pcobip, Mode="s")
pcobip=AddCluster2Biplot(pcobip, NGroups=3, ClusterType="hi")
op <- par(mfrow=c(1,2)) 
plot(pcobip, Mode="s", PlotClus = TRUE)
plot(pcobip$Dendrogram)
par(op)

}