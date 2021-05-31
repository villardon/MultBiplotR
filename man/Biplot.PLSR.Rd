\name{Biplot.PLSR}
\alias{Biplot.PLSR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Biplot
}
\description{
Adds a Biplot to a Partial Lest Squares (plsr) object.
}
\usage{
Biplot.PLSR(plsr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plsr}{
A plsr object from the PLSR function
}
}
\details{
Adds a Biplot to a Partial Lest Squares (plsr) object. The biplot is constructed with the matrix of predictors, the dependent variable is projected onto the biplot as a continuous supplementary variable.
}
\value{
An object of class ContinuousBiplot with the dependent variables as supplemntary.
}
\references{
Oyedele, O. F., & Lubbe, S. (2015). The construction of a partial least-squares biplot. Journal of Applied Statistics, 42(11), 2449-2460.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{PLSR}}
}
\examples{
\donttest{X=as.matrix(wine[,4:21])
y=as.numeric(wine[,2])-1
mifit=PLSR(y,X, Validation="None")
mibip=Biplot.PLSR(mifit)
plot(mibip, PlotVars=TRUE, IndLabels = y, ColorInd=y+1)}
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
