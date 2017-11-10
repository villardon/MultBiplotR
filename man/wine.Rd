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
    \item{\code{Year}}{A factor with levels \code{1986} \code{1987}}
    \item{\code{Origin}}{A factor with levels \code{Ribera} \code{Toro}}
    \item{\code{Group}}{A factor with levels \code{R86} \code{R87} \code{T86} \code{T87}}
    \item{\code{A}}{Alcoholic content (percentage)}
    \item{\code{VA}}{volatil acidity - g acetic acid/l}
    \item{\code{TA}}{Total tritable acidity - g tartaric acid/l}
    \item{\code{FA}}{Fixed acidity -  g tartaric acid/l}
    \item{\code{pH}}{ph}
    \item{\code{TPR}}{Total phenolics - g gallic acid /l - Folin}
    \item{\code{TPS}}{Total phenolics - Somers }
    \item{\code{V}}{Substances reactive to vanilin - mg catechin/l}
    \item{\code{PC}}{Procyanidins - mg cyanidin/l}
    \item{\code{ACR}}{Total Anthocyanins - mg/l - method 1}
    \item{\code{ACS}}{Total Anthocyanins - mg/l - methods 2}
    \item{\code{ACC}}{Malvidin - malvidin-3-glucoside mg/l}
    \item{\code{CI}}{Color density - }
    \item{\code{CI2}}{Color density 2}
    \item{\code{H}}{Wine Hue Color}
    \item{\code{I}}{Degree of Ionization - Percent}
    \item{\code{CA}}{Chemical Age}
    \item{\code{VPC}}{ratio V/PC}
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
