\name{CanonicalStatisBiplot}
\alias{CanonicalStatisBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CANONICAL STATIS-ACT for multiple tables with common rows and its associated Biplot
}
\description{
The procedure performs STATIS-ACT methodology for multiple tables with common rows and its associated biplot
}
\usage{
CanonicalStatisBiplot(X, Groups, InitTransform = "Standardize columns", dimens = 2,
                 SameVar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A list containing multiple tables with common rows
}
  \item{Groups}{
A factor containing the groups
}
  \item{InitTransform}{
Initial transformation of the data matrices
}
  \item{dimens}{
Dimension of the final solution
}
  \item{SameVar}{
Are the variables the same for all occasions?
}
}
\details{
The procedure performs Canonical STATIS-ACT methodology for multiple tables with common rows and its associated biplot. When the variables are the same for all occasions trajectories for the variables can also be plotted.
}
\value{
An object of class StatisBiplot
}
\references{
Vallejo-Arboleda, A., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2007). Canonical STATIS: Biplot analysis of multi-table group structured data based on STATIS-ACT methodology. Computational statistics & data analysis, 51(9), 4193-4205.

Abdi, H., Williams, L.J., Valentin, D., & Bennani-Dosse, M. (2012). STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat, 4, 124-167.

Efron, B.,Tibshirani, RJ. (1993). An introduction to the bootstrap. New York: Chapman and Hall. 436p.

Escoufier, Y. (1976). Operateur associe a un tableau de donnees. Annales de laInsee, 22-23, 165-178.

Escoufier, Y. (1987). The duality diagram: a means for better practical applications. En P. Legendre & L. Legendre (Eds.), Developments in Numerical Ecology, pp. 139-156, NATO Advanced Institute, Serie G. Berlin: Springer.

L'Hermier des Plantes, H. (1976). Structuration des Tableaux a Trois Indices de la Statistique. [These de Troisieme Cycle]. University of Montpellier, France.

Ringrose, T.J. (1992). Bootstrapping and Correspondence Analysis in Archaeology. Journal of Archaeological Science. 19:615-629.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{data(Chemical)
x= Chemical[37:144,5:9]
weeks=as.factor(as.numeric(Chemical$WEEKS[37:144]))
levels(weeks)=c("W2" , "W3", "W4")
X=Convert2ThreeWay(x,weeks, columns=FALSE)
Groups=Chemical$Treatment[1:36]
canstbip=CanonicalStatisBiplot(X, Groups, SameVar = TRUE)
plot(canstbip, mode="s", PlotVars=TRUE, ShowBox=TRUE)}
}

