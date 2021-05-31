\name{CorrelationCircle}
\alias{CorrelationCircle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Circle of correlations
}
\description{
Circle of correlations among the manifiest variables and the principal comonents (or dimensions of any biplot).
}
\usage{
CorrelationCircle(bip, A1 = 1, A2 = 2, Colors = NULL, Labs = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
an biplot object of any kind.
}
  \item{A1}{
First dimension for the representation
}
  \item{A2}{
Second dimension for the representation
}
  \item{Colors}{
Colors of the variables
}
  \item{Labs}{
Labels of the variables
}
  \item{\dots}{
Any other graphical parameters
}
}
\details{
Circle of correlations among the manifiest variables and the principal comonents (or dimensions of any biplot).
}
\value{
The plot of the circle of correlations
}

\author{
Jose Luis Vicente Villardon
}

\examples{
bip=PCA.Biplot(wine[,4:21])
CorrelationCircle(bip)
}
