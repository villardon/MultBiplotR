\name{AddContVars2Biplot}
\alias{AddContVars2Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds supplementary continuous variables to a biplot object
}
\description{
Adds supplementary continuous variables to a biplot object
}
\usage{
AddContVars2Biplot(bip, X, dims = NULL, Scaling = 5, Fit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
A biplot object
}
  \item{X}{
Matrix containing the supplementary continuos variables
}
  \item{dims}{
Dimension of the solution
}
  \item{Scaling}{
Transformation to apply to X
}
  \item{Fit}{
Type of fit. Linear by default.
}
}
\details{
More types of fit will be added in the future
}
\value{
A biplot object with the  coordinates for the supplementary variables added.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{AddSupVars2Biplot}}
}
\examples{
# Not yet

}