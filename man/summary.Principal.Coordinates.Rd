\name{summary.Principal.Coordinates}
\alias{summary.Principal.Coordinates}
\title{Summary of the results of a Principal Coordinates Analysis}


\usage{
\method{summary}{Principal.Coordinates}(object, printdata=FALSE, printproximities=FALSE, 
printcoordinates=FALSE, printqualities=FALSE,\dots)
}

\arguments{
  \item{object}{
An object of Type \code{Principal.Coordinates}
}
  \item{printdata}{
Should original data be printed. Default is FALSE
}
  \item{printproximities}{
Should proximities be printed. Default is FALSE
}
  \item{printcoordinates}{
Should proximities be printed. Default is FALSE
}
  \item{printqualities}{
Should qualoties of representation be printed. Default is FALSE
}
  \item{\dots}{
Additional parameters to summary.
}
}

\description{
  Summary of the results of a Principal Coordinates Analysis
}


\details{
This function is a method for the generic function summary() for class "Principal.Coordinates". It can be invoked by calling summary(x) for an object x of the appropriate class.
}

\value{
The summary
}
\author{Jose Luis Vicente-Villardon
}

\examples{
data(spiders)
dist=BinaryProximities(spiders)
pco=PrincipalCoordinates(dist)
summary(pco)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{summary}
