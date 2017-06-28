\name{plot.CA.sol}
\alias{plot.CA.sol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the solution of a Coorespondence Analysis
}
\description{
Plots the solution of a Correspondence Analysis
}
\usage{
\method{plot}{CA.sol}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A CA.sol object
}
  \item{\dots}{
Any other biplot and graphical parameters
}
}
\details{
Plots the solution of a Correspondence Analysis
}
\value{
No value returned
}
\references{
Add some references here
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
data(riano)
Sp=riano[,3:15]
cabip=CA(Sp)
plot(cabip)
}

