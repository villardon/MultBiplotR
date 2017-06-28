\name{plot.Statis}
\alias{plot.Statis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a Statis Object
}
\description{
Plots a Statis Object
}
\usage{
\method{plot}{Statis}(x, A1 = 1, A2 = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A Statis object
}
  \item{A1}{
First dimension of the plot
}
  \item{A2}{
Second dimension of the plot
}
  \item{\dots}{
Aditional parameters}
}
\details{
Plots a Statis Object
}
\value{
A biplot
}
\references{
Vallejo-Arboleda, A., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2007). Canonical STATIS: Biplot analysis of multi-table group structured data based on STATIS-ACT methodology. Computational statistics & data analysis, 51(9), 4193-4205.
}
\author{
Jose Luis Vicente Villardon
}
\examples{
data(Chemical)
x= Chemical[,5:16]
X=Convert2ThreeWay(x,Chemical$WEEKS, columns=FALSE)
stbip=StatisBiplot(X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
