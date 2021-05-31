\name{CCA}
\alias{CCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Canonical Correspondence Analysis
}
\description{
Calculates the solution of a Canonical Correspondence Analysis Biplot
}
\usage{
CCA(P, Z, alpha = 1, dimens = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{P}{
Abundance Matrix of sites by species.
}
  \item{Z}{
Environmental variables measured at the same sites
}
  \item{alpha}{
Alpha for the biplot decomposition [0,1]. With alpha=1 the emphasis is on the sites
    and with alpha=0 the emphasis is on the species
}
  \item{dimens}{
Dimension of the solution
}

}
\details{
A pair of ecological tables, made of a species abundance matrix and an environmental variables matrix measured at the same sampling sites, is usually analyzed by Canonical Correspondence Analysis (CCA) (Ter BRAAK, 1986).  CCA can be considered as a Correspondence Analysis (CA) in which the ordination axis are constrained to be linear combinations of the environmental variables. Recently the procedure has been extended to other disciplines as Sociology or Psichology and it is potentially useful in many other fields.
}
\value{
A CCA solution object
}
\references{
Ter Braak, C. J. (1986). Canonical correspondence analysis: a new eigenvector technique for multivariate direct gradient analysis. Ecology, 67(5), 1167-1179.

Johnson, K. W., & Altman, N. S. (1999). Canonical correspondence analysis as an approximation to Gaussian ordination. Environmetrics, 10(1), 39-52.

Graffelman, J. (2001). Quality statistics in canonical correspondence analysis. Environmetrics, 12(5), 485-497.

Graffelman, J., & Tuft, R. (2004). Site scores and conditional biplots in canonical correspondence analysis. Environmetrics, 15(1), 67-80.

Greenacre, M. (2010). Canonical correspondence analysis in social science research (pp. 279-286). Springer Berlin Heidelberg.

}
\author{
Jose Luis vicente Villardon
}

\examples{
\donttest{data(riano)
Sp=riano[,3:15]
Env=riano[,16:25]
ccabip=CCA(Sp, Env)
plot(ccabip)}
}

