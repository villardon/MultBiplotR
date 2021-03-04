\name{StatisBiplot}
\alias{StatisBiplot}

\title{
STATIS-ACT for multiple tables with common rows and its associated Biplot
}
\description{
The procedure performs STATIS-ACT methodology for multiple tables with common rows and its associated biplot
}
\usage{
StatisBiplot(X, InitTransform = "Standardize columns", dimens = 2,
                 SameVar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A list containing multiple tables with common rows.
}
  \item{InitTransform}{
Initial transformation of the data matrices
}
  \item{dimens}{
Dimension of the final solution
}
  \item{SameVar}{
Are the variables the same for all occasions? If so, Biplot trajectories for each variable will be calculated.
}
}
\details{
The procedure performs STATIS-ACT methodology for multiple tables with common rows and its associated biplot. When the variables are the same for all occasions trajectories for the variables can also be plotted.
Basic plotting includes the consensus individuals and all the variables. Traditional trajectories for individuals and biplot trajectories for variables (when adequate) are optional.
The original matrix will be provided as a list each cell of the list is the data matrix for one ocassion
the number of rows for each occasion must be the same
}
\value{
An object of class StatisBiplot
}
\references{
Abdi, H., Williams, L.J., Valentin, D., & Bennani-Dosse, M. (2012). STATIS and DISTATIS: optimum multitable principal component analysis and three way metric multidimensional scaling. WIREs Comput Stat, 4, 124-167.

Efron, B.,Tibshirani, RJ. (1993). An introduction to the bootstrap. New York: Chapman and Hall. 436p.

Escoufier, Y. (1976). Operateur associe a un tableau de donnees. Annales de laInsee, 22-23, 165-178.

Escoufier, Y. (1987). The duality diagram: a means for better practical applications. En P. Legendre & L. Legendre (Eds.), Developments in Numerical Ecology, pp. 139-156, NATO Advanced Institute, Serie G. Berlin: Springer.

L'Hermier des Plantes, H. (1976). Structuration des Tableaux a Trois Indices de la Statistique. [These de Troisieme Cycle]. University of Montpellier, France.

Ringrose, T.J. (1992). Bootstrapping and Correspondence Analysis in Archaeology. Journal of Archaeological. Science.19:615-629.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(Chemical)
# Extract continous data from the original data frame.
x= Chemical[,5:16]
# Obtaining the three way table as a list
X=Convert2ThreeWay(x,Chemical$WEEKS, columns=FALSE)
# Calculating the Biplot associated to STATIS-ACT
stbip=StatisBiplot(X, SameVar=TRUE)
# Basic plot of the results
plot(stbip)
# Colors By Table
plot(stbip, VarColorType="ByTable")
# Colors By Variable
plot(stbip, VarColorType="ByVar", mode="s",  MinQualityVars = 0.5)

plot(stbip, PlotRowTraj = TRUE, PlotVars=FALSE, RowColors=1:36)
}

