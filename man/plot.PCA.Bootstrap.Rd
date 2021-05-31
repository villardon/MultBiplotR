\name{plot.PCA.Bootstrap}
\alias{plot.PCA.Bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the Bootstrap information for Principal Components Analysis (PCA)
}
\description{
Plots an object of class "PCA.Bootstrap"
}

\usage{
\method{plot}{PCA.Bootstrap}(x, Eigenvalues = TRUE, 
Inertia = FALSE, EigenVectors = TRUE, Structure = TRUE, 
Squared = TRUE, Scores = TRUE, ColorInd = "black", TypeScores = "ch", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class "PCA.Bootstrap"
}
  \item{Eigenvalues}{
Should the information for the eigenvalues be plotted?
}
  \item{Inertia}{
Should the information for the inertia be plotted?
}
  \item{EigenVectors}{
Should the information for the eigenvectors be plotted?
}
  \item{Structure}{
Should the information for the correlations (variables-dimensions) be plotted?
}
  \item{Squared}{
Should the information for the correlations (variables-dimensions) be plotted?
}
  \item{Scores}{
Should the row (individual) scores be plotted?
}
  \item{ColorInd}{
Colors for the rows
}
  \item{TypeScores}{
Type of plot for the scores
}
  \item{...}{
Any other graphical argument
}
}
\details{
For each parameter, box-plots and confidence intervals are plotted. 
The initial estimator and the bootstrap mean are plotted.

For the eigenvectors, loadings and contributions, the graph is divided into as many rows as dimensions, each row contains a plot of the hole set of variables.

The scores are plotted on a two dimensional


}
\value{
No value returned 
}
\references{

Daudin, J. J., Duby, C., & Trecourt, P. (1988). Stability of principal component analysis studied by the bootstrap method. Statistics: A journal of theoretical and applied statistics, 19(2), 241-258.

Chateau, F., & Lebart, L. (1996). Assessing sample variability in the visualization techniques related to principal component analysis: bootstrap and alternative simulation methods. COMPSTAT, Physica-Verlag, 205-210.

Babamoradi, H., van den Berg, F., & Rinnan, Å. (2013). Bootstrap based confidence limits in principal component analysis: A case study. Chemometrics and Intelligent Laboratory Systems, 120, 97-105.

Fisher, A., Caffo, B., Schwartz, B., & Zipunnikov, V. (2016). Fast, exact bootstrap principal component analysis for p> 1 million. Journal of the American Statistical Association, 111(514), 846-860.


}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{PCA.Bootstrap}}
}
\examples{
\donttest{X=wine[,4:21]
grupo=wine$Group
rownames(X)=paste(1:45, grupo, sep="-")
pcaboot=PCA.Bootstrap(X, dimens=2, Scaling = "Standardize columns", B=1000)
plot(pcaboot, ColorInd=as.numeric(grupo))
summary(pcaboot)}
}

