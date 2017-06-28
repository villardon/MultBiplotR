\name{wine}
\alias{wine}
\docType{data}
\title{
Wine data
}
\description{
Comparison of young wines of Ribera de Duero and Toro
}
\usage{data("wine")}
\format{
  A data frame with 45 observations on the following 21 variables.
  \describe{
    \item{\code{Year}}{a factor with levels \code{1986} \code{1987}}
    \item{\code{Origin}}{a factor with levels \code{Ribera} \code{Toro}}
    \item{\code{Group}}{a factor with levels \code{R86} \code{R87} \code{T86} \code{T87}}
    \item{\code{A}}{a numeric vector}
    \item{\code{VA}}{a numeric vector}
    \item{\code{TA}}{a numeric vector}
    \item{\code{FA}}{a numeric vector}
    \item{\code{pH}}{a numeric vector}
    \item{\code{TPR}}{a numeric vector}
    \item{\code{TPS}}{a numeric vector}
    \item{\code{V}}{a numeric vector}
    \item{\code{PC}}{a numeric vector}
    \item{\code{ACR}}{a numeric vector}
    \item{\code{ACS}}{a numeric vector}
    \item{\code{ACC}}{a numeric vector}
    \item{\code{CI}}{a numeric vector}
    \item{\code{CI2}}{a numeric vector}
    \item{\code{H}}{a numeric vector}
    \item{\code{I}}{a numeric vector}
    \item{\code{CA}}{a numeric vector}
    \item{\code{VPC}}{a numeric vector}
  }
}
\details{
Comparison of young wines of Ribera de Duero and Toro
}
\source{
Rivas-Gonzalo, J. C., Gutierrez, Y., Polanco, A. M., Hebrero, E., Vicente-Villardon, J. L., Galindo, P., & Santos-Buelga, C. (1993). Biplot analysis applied to enological parameters in the geographical classification of young red wines. American journal of enology and viticulture, 44(3), 302-308.
}
\references{
Rivas-Gonzalo, J. C., Gutierrez, Y., Polanco, A. M., Hebrero, E., Vicente-Villardon, J. L., Galindo, P., & Santos-Buelga, C. (1993). Biplot analysis applied to enological parameters in the geographical classification of young red wines. American journal of enology and viticulture, 44(3), 302-308.
}
\examples{
data(wine)
## maybe str(wine) ; plot(wine) ...
}
\keyword{datasets}
