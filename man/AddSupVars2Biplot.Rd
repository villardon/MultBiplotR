\name{AddSupVars2Biplot}
\alias{AddSupVars2Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds supplementary variables to a biplot object
}
\description{
Adds supplementary bariables to a biplot object constructed with any of the biplot methods of the package. The new variables are fitted using the coordinates for the rows. Each variable is fitted using the adequate procedure for its type.
}
\usage{
AddSupVars2Biplot(bip, X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
The biplot object
}
  \item{X}{
A data frame with the supplementary variables.
}
}
\details{
Binary, nominal or ordinal variables are fitted using logistic biplots. Continuous variables are fitted with linear regression.
}
\value{
A biplot object with the  coordinates for the supplementary variables added.
}

\author{
Jose Luis Vicente Villardon
}

\seealso{
 \code{\link{AddContVars2Biplot}}
}
\examples{
# Not yet

}