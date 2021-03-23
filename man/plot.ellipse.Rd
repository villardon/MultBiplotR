\name{plot.ellipse}
\alias{plot.ellipse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a concentration ellipse.
}
\description{
Plot a concentration ellipse obtained from \code{ConcEllipse}.
}
\usage{
\method{plot}{ellipse}(x, add=TRUE, labeled= FALSE , 
center=FALSE, centerlabel="Center", initial=FALSE,  ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object with class \code{ellipse} obtained from \code{ConcEllipse}.
}
  \item{add}{
Should the ellipse be added to the current plot?
}
  \item{labeled}{
Should the ellipse be labelled with the confidence level?
}
  \item{center}{
Should the center be plotted?
}

  \item{centerlabel}{
Label for the center.
}

  \item{initial}{
Should the initial data be plotted?
}

  \item{\dots}{
Any other graphical parameter that can affects the plot (as color, etc ...)
}
}
\details{
Plots an ellipse containing a specified percentage of the data.
}
\value{
No value returned
}
\references{
Meulman, J. J., & Heiser, W. J. (1983). The display of bootstrap solutions in multidimensional scaling. Murray Hill, NJ: Bell Laboratories.

Linting, M., Meulman, J. J., Groenen, P. J., & Van der Kooij, A. J. (2007). Stability of nonlinear principal components analysis: An empirical study using the balanced bootstrap. Psychological Methods, 12(3), 359.
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{ConcEllipse}}, ~~~
}
\examples{
data(iris)
dat=as.matrix(iris[1:50,1:2])
plot(iris[,1], iris[,2],col=iris[,5], asp=1)
E=ConcEllipse(dat, 0.95)
plot(E, labeled=TRUE, center=TRUE)
}

