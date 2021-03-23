\name{plot.fraction}
\alias{plot.fraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a fraction of the data as a cluster
}
\description{
Plots a convex hull or a star  containing a specified percentage of the data. Used to plot clusters.
}
\usage{
\method{plot}{fraction}(x, add = TRUE, center = FALSE, 
centerlabel = "Center", initial = FALSE, type = "ch", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object with class \code{fraction} obtained from \code{Fraction}.
}
  \item{add}{
Should the fraction be added to the current plot?
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
  \item{type}{
Type of plot. Can be: "ch"- Convex Hull or "st" - Star (Joining each point with the center)
}
  \item{\dots}{
Any other graphical parameter that can affects the plot (as color, etc ...)
}
}
\details{
Plots a convex hull or a star  containing a specified percentage of the data.
}
\value{
No value returned
}

\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Fraction}}
}
\examples{
a=matrix(runif(50), 25,2)
a2=Fraction(a, 0.7)
plot(a2, add=FALSE, type="ch", initial=TRUE, center=TRUE, col="blue")
plot(a2, add=TRUE, type="st", col="red")
}

