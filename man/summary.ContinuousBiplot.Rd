\name{summary.ContinuousBiplot}
\alias{summary.ContinuousBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary of the solution of a Biplot for Continuous Data
}
\description{
Summary of the solution of a Biplot for Continuous Data
}
\usage{
\method{summary}{ContinuousBiplot}(object, latex = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class "ContinuousBiplot"
}
  \item{latex}{
Should the results be in latex tables
}
  \item{\dots}{
Any aditional parameters
}
}
\details{
Summary of the solution of a Biplot for Continuous Data
}
\value{
The summary
}
\author{
Jose Luis Vicente Villardon
}
\examples{
## Simple Biplot with arrows
data(Protein)
bip=PCA.Biplot(Protein[,3:11])
summary(bip)
}

