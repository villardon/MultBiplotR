\name{AddOrdVars2Biplot}
\alias{AddOrdVars2Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adds supplementary ordinal variables to an existing biplot objects.
}
\description{
Adds supplementary ordinal variables to an existing biplot objects.
}
\usage{
AddOrdVars2Biplot(bip, Y, tol = 1e-06, maxiterlogist = 100, 
penalization = 0.2, showiter = TRUE, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
A biplot object.
}
  \item{Y}{
A matrix of ordinal variables.
}
  \item{tol}{
Tolerance.
}
  \item{maxiterlogist}{
Maximum number of iterations for the logistic fit.
}
  \item{penalization}{
Penalization for the logistic fit
}
  \item{showiter}{
Should the itrations be shown on screen
}
  \item{show}{
Show details.
}
}
\details{
Adds supplementary ordinal variables to an existing biplot objects.
}
\value{
An object with the information of the fits
}
\references{
Vicente-Villardon, J. L., & Hernandez-Sanchez, J. C. (2020). External Logistic Biplots for Mixed Types of Data. In Advanced Studies in Classification and Data Science (pp. 169-183). Springer, Singapore.
}
\author{
Jose Luis Vicente-Villardon
}

\examples{
# not yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
