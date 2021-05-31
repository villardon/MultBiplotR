\name{AddCluster2Biplot}
\alias{AddCluster2Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Add clusters to a biplot object
}
\description{
The function add clusters to a biplot object to be represented on the biplot. The clusters can be defined by a nominal variable provided by the user, obtained from the \code{hclust} function of the base package or from the \code{kmeans} function
}
\usage{
AddCluster2Biplot(Bip, NGroups=3, ClusterType="hi", Groups=NULL, Original=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bip}{
A Biplot object obtained from any biplot procedure. It has to be a list containing a field called \code{Bip$RowCoordinates} in order to calculate the clusters when necessary. 
}
  \item{NGroups}{
Number of groups or clusters. Only necessary when hierarchical or k-means procedures are used.
}
  \item{ClusterType}{
The type of cluster to add. There are three possibilities "us" (User Defined), "hi" (hierarchical clusters),  "km" (kmeans clustering) or "gm" (gaussian mixture).
}
  \item{Groups}{
A factor defining the groups provided by the user.
}
  \item{Original}{
Should the clusters be calculated using the original data rather than the reduced dimensions?.
}
  \item{...}{
Any other parameter for the \code{hclust} and \code{kmeans} procedures.
}
}
\details{
One of the main shortcomings of cluster analysis is that it is not easy to search for the variables associated to the obtained classification; representing the clusters on the biplot can help to perform that interpretation. If you consider the technique for dimension reduction as a way to separate the signal from the noise, clusters should be constructed using the dimensions retained in the biplot, otherwise the complete original data matrix can be used. The colors used by each cluster should match the color used in the Dendrogram.
User defined clusters can also be plotted, for example, to investigate the relation of the biplot solution to an external nominal variable.
}
\value{
The function returns the biplot object with the information about the clusters added in new fields
\item{ClusterType }{The method of clustering as defined in the argument \code{ClusterType}.}
\item{Clusters}{A factor containing the solution or the user defined clusters}
\item{ClusterNames}{The names of the clusters}
\item{ClusterColors}{The colors of the clusters}
\item{Dendrogram}{The Dendrogram if we have used hirarchical clustering}
\item{ClusterObject}{The object obtained from \code{hclust}, \code{kmeans} or \code{MGC}}
}
\references{
Demey, J. R., Vicente-Villardon, J. L., Galindo-Villardon, M. P., & Zambrano, A. Y. (2008). Identifying molecular markers associated with classification of genotypes by External Logistic Biplots. Bioinformatics, 24(24), 2832-2838.

Gallego-Alvarez, I., & Vicente-Villardon, J. L. (2012). Analysis of environmental indicators in international companies by applying the logistic biplot. Ecological Indicators, 23, 250-261.

Galindo, P. V., Vaz, T. D. N., & Nijkamp, P. (2011). Institutional capacity to dynamically innovate: an application to the Portuguese case. Technological Forecasting and Social Change, 78(1), 3-12.

Vazquez-de-Aldana, B. R., Garcia-Criado, B., Vicente-Tavera, S., & Zabalgogeazcoa, I. (2013). Fungal Endophyte (Epichloë festucae) Alters the Nutrient Content of Festuca rubra Regardless of Water Availability. PloS one, 8(12), e84539.

}
\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
For clusters not provided by the user the function uses the standard procedures in  \code{hclust} and \code{kmeans}.
}
\examples{
\donttest{
data(Protein)
bip=PCA.Biplot(Protein[,3:11])
plot(bip)
# Add user defined clusters containing the region (North, South, Center)
bip=AddCluster2Biplot(bip, ClusterType="us", Groups=Protein$Region)
plot(bip, mode="a", margin=0.1, PlotClus=TRUE)

# Hierarchical clustering on the biplot coordinates using the Ward method
bip=AddCluster2Biplot(bip, ClusterType="hi", method="ward.D")
op <- par(mfrow=c(1,2))
plot(bip, mode="s", margin=0.1, PlotClus=TRUE)
plot(bip$Dendrogram)
par(op)
# K-means cluster on the biplot coordinates using the Ward method
bip=AddCluster2Biplot(bip, ClusterType="hi", method="ward.D")
op <- par(mfrow=c(1,2))
plot(bip, mode="s", margin=0.1, PlotClus=TRUE)
plot(bip$Dendrogram)
par(op)
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Cluster Analysis }
\keyword{ Biplot }% __ONLY ONE__ keyword per line
