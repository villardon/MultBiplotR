\name{OrdLogBipEM}
\alias{OrdLogBipEM}
\title{
Alternated EM algorithm for Ordinal Logistic Biplots
}
\description{
This function computes, with an alternated algorithm, the row and column parameters of an Ordinal Logistic Biplot for ordered polytomous data. The row coordinates (E-step) are computed using multidimensional Gauss-Hermite quadratures and Expected \emph{a posteriori} (EAP) scores and parameters for each variable or items (M-step) using Ridge Ordinal Logistic Regression to solve the separation problem present when the points for different categories of a variable are completely separated on the representation plane and the usual fitting methods do not converge. The separation problem is present in almost avery data set for which the goodness of fit is high. 
}
          
\usage{
OrdLogBipEM(Data, freq=NULL, dim = 2, nnodes = 15, 
tol = 0.0001, maxiter = 100, maxiterlogist = 100, 
penalization = 0.2, show = FALSE, initial = 1, alfa = 1, 
Orthogonalize=TRUE, Varimax=TRUE, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
Data frame with the ordinal data. All the variables must be ordered factors. 
}
  \item{freq}{
Frequencies for compacted tables
}
  \item{dim}{
Dimension of the solution
}
  \item{nnodes}{
Number of nodes for the multidimensional Gauss-Hermite quadrature
}
 \item{tol}{
Value to stop the process of iterations.
}
  \item{maxiter}{
Maximum number of iterations for the biplot procedure.
}
  \item{maxiterlogist}{
Maximum number of iterations for the logistic regression step or the Mirt initial configuration.
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
\item{show}{
Boolean parameter to specify if the user wants to see every iteration.
}
\item{initial}{
Method used to choose the initial ability in the algorithm. Default value is 1.
}
\item{alfa}{
Optional parameter to calculate row and column coordinates in Simple correspondence analysis if the initial parameter is equal to 1.
}
\item{Orthogonalize}{
Should the final row coordinates be orthogonalized?. The column parameters have to be recalculated.
}
\item{Varimax}{
Should the final row coordinates be rotated using the varimax procedure?.
}
\item{...}{
Aditional argunments for mirt.
}
}

\value{
An object of class \code{"Ordinal.Logistic.Biplot"}.This has components:

  \item{RowCoordinates}{Coordinates for the rows or the individuals}
  \item{ColumnParameters}{List with information about the Ordinal Logistic Models calculated for each variable including: estimated parameters with thresholds,percents of correct classifications,and pseudo-Rsquared}
  \item{loadings}{factor loadings}
  \item{LogLikelihood}{Logarithm of the likelihood}
  \item{r2}{R squared coefficient}
  \item{Ncats}{Number of the categories of each variable}

}
\references{
  Bock,R. & Aitkin,M. (1981),\emph{Marginal maximum likelihood estimation of item parameters: Aplication of an EM algorithm}, Phychometrika 46(4), 443-459.
  }
  
\author{
Jose Luis Vicente-Villardon
}

\examples{
\dontrun{
    data(Doctors)
    olb = OrdLogBipEM(Doctors,dim = 2, nnodes = 10, initial=4,
    tol = 0.001, maxiter = 100, penalization = 0.1, show=TRUE)
    olb
    summary(olb)
    PlotOrdinalResponses(olb)
    }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EM}
\keyword{algorithm}