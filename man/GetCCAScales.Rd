\name{GetCCAScales}
\alias{GetCCAScales}

\title{
Calculates scales for plotting the environmental variables in a Canonical Correspondence Analysis
}
\description{
Calculates scales for plotting the environmental variables in a Canonical Correspondence Analysis
}
\usage{
GetCCAScales(CCA, nticks = 7, TypeScale = "Complete", ValuesScale = "Original")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CCA}{
A CCA solution object
}
  \item{nticks}{
Number of ticks to represent
}
  \item{TypeScale}{
Type of scale to represent
}
  \item{ValuesScale}{
Values to represent (Original or Transformed)
}
}
\details{
Calculates scales for plotting the environmental variables in a Canonical Correspondence Analysis
}
\value{
Returns the points and the labels for each biplot axis
}
\references{
Gower, J. C., & Hand, D. J. (1995). Biplots (Vol. 54). CRC Press.

Gower, J. C., Lubbe, S. G., & Le Roux, N. J. (2011). Understanding biplots. John Wiley & Sons.

Vicente-Villardón, J. L., Galindo Villardón, M. P., & Blázquez Zaballos, A. (2006). Logistic biplots. Multiple correspondence analysis and related methods. London: Chapman & Hall, 503-521.

}
\author{
Jose Luis Vicente Villardon
}

\examples{
# No examples yet
}
