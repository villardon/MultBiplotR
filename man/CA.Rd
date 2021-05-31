\name{CA}
\alias{CA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correspondence Analysis
}
\description{
Correspondence Analysis for a frequency or abundace data matrix.
}
\usage{
CA(x, dim = 2, alpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The frequency or abundance data matrix.
}
  \item{dim}{
Dimension of the final solution
}
  \item{alpha}{
Alpha to determine the kind of biplot to use.
}
}
\details{
Calculates Correspondence Analysis for a tww-way frequency or abundance table
}
\value{
Correspondence analysis solution
}
\references{
Benzécri, J. P. (1992). Correspondence analysis handbook. New York: Marcel Dekker.

Greenacre, M. J. (1984). Theory and applications of correspondence analysis. Academic Press.
}
\author{
Jose Luis Vicente Villardon
}


\examples{
\donttest{data(SpidersSp)
cabip=CA(SpidersSp)
plot(cabip)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

