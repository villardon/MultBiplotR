\name{ContinuousDistances}
\alias{ContinuousDistances}
\title{Distances for Continuous Data}
\description{
Calculates distances among rows of a continuous data matrix or among the rows of two continuous matrices.
}
\usage{
  ContinuousDistances(x, y = NULL, coef = "Pythagorean", r = 1)
}

\arguments{
  \item{x}{
Main data matrix. Distances among rows are calculated if y=NULL.
}
  \item{y}{
Supplementary data matrix. If not NULL the distances among the rows of x and y are calculated 
}
  \item{coef}{
Distance coefficient. Use the name or the number(see details)
}
  \item{r}{
Exponent for the Minkowsky
}
}
\details{
The following coefficients are calculated

1.- Pythagorean =  sqrt(sum((y[i, ] - x[j, ])^2)/p)

2.- Taxonomic =  sqrt(sum(((y[i,]-x[j,])^2)/r^2)/p)

3.- City = sum(abs(y[i,]-x[j,])/r)/p

4.- Minkowski = (sum((abs(y[i,]-x[j,])/r)^t)/p)^(1/t)

5.- Divergence = sqrt(sum((y[i,]-x[j,])^2/(y[i,]+x[j,])^2)/p)

6.- dif_sum = sum(abs(y[i,]-x[j,])/abs(y[i,]+x[j,]))/p

7.- Camberra = sum(abs(y[i,]-x[j,])/(abs(y[i,])+abs(x[j,])))

8.- Bray_Curtis = sum(abs(y[i,]-x[j,]))/sum(y[i,]+x[j,])

9.- Soergel = sum(abs(y[i,]-x[j,]))/sum(apply(rbind(y[i,],x[j,]),2,max))

10.- Ware_hedges = sum(abs(y[i,]-x[j,]))/sum(apply(rbind(y[i,],x[j,]),2,max))
}

\value{A list with:
\item{Data}{A matrix with the initial data (x matrix).}
\item{SupData}{A matrix with the supplementary data (y matrix).}
\item{D}{The matrix of distances}
\item{Coefficient}{The coefficient used.}
}

\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley
}

\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{PrincipalCoordinates}}
}
\examples{
\donttest{data(wine)
dis=ContinuousDistances(wine[,4:21])}
}
