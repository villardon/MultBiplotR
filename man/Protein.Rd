\name{Protein}
\alias{Protein}
\docType{data}
\title{
Protein consumption data.
}
\description{
Protein consumption in twenty-five European countries for nine food groups.
}
\usage{data(Protein)}
\format{
  A data frame with 25 observations on the following 11 variables.
  \describe{
    \item{\code{Comunist}}{a factor with levels \code{No} \code{Yes}}
    \item{\code{Region}}{a factor with levels \code{North} \code{Center} \code{South}}
    \item{\code{Red_Meat}}{a numeric vector}
    \item{\code{White_Meat}}{a numeric vector}
    \item{\code{Eggs}}{a numeric vector}
    \item{\code{Milk}}{a numeric vector}
    \item{\code{Fish}}{a numeric vector}
    \item{\code{Cereal}}{a numeric vector}
    \item{\code{Starch}}{a numeric vector}
    \item{\code{Nuts}}{a numeric vector}
    \item{\code{Fruits_Vegetables}}{a numeric vector}
  }
}
\details{
These data measure protein consumption in twenty-five European countries for nine food groups. It is possible to use multivariate methods to determine whether there are groupings of countries and whether meat consumption is related to that of other foods.
}
\source{
http://lib.stat.cmu.edu/DASL/Datafiles/Protein.html
}
\references{
Weber, A. (1973) Agrarpolitik im Spannungsfeld der internationalen Ernaehrungspolitik, Institut fuer Agrarpolitik und marktlehre, Kiel. 

Gabriel, K.R. (1981) Biplot display of multivariate matrices for inspection of data and diagnosis. In Interpreting Multivariate Data (Ed. V. Barnett), New York: John Wiley & Sons, 147-173. 

Hand, D.J., et al. (1994) A Handbook of Small Data Sets, London: Chapman & Hall, 297-298.
}
\examples{
data(Protein)
## maybe str(Protein) ; plot(Protein) ...
}
\keyword{datasets}
