\name{OrdVarBiplot}
\alias{OrdVarBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots an ordinal variable on the biplot
}
\description{
Plots an ordinal variable on the biplot from its fitted parameters
}
\usage{
OrdVarBiplot(bi1, bi2, threshold, xmin = -3, xmax = 3, ymin = -3, 
ymax = 3, label = "Point", mode = "a", CexMarks = 0.7, CexPoint = 0.8, 
PchPoint = 1, Color = "green", tl = 0.03, textpos = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bi1}{
Slope for the first dimension to plot
}
  \item{bi2}{
Slope for the second dimension to plot
}
  \item{threshold}{
Thresholds for each category of the variable}
  \item{xmin}{
Minimum value of the X on the plot
}
  \item{xmax}{
Maximum value of the X on the plot
}
  \item{ymin}{
Minimum value of the Y on the plot
}
  \item{ymax}{
Maximum value of the X on the plot
}
  \item{label}{
Label of the variable
}
  \item{mode}{
Mode of the plot (as in a regular biplot)
}
  \item{CexMarks}{
Size of the tick marks
}
  \item{CexPoint}{
Size of the point
}
  \item{PchPoint}{
Mark for the point
}
  \item{Color}{
Color
}
  \item{tl}{
Tick Length
}
  \item{textpos}{
Position of the label
}
  \item{\dots}{
Any aditional graphical parameter
}
}
\details{
Plots an ordinal variable on the biplot from its fitted parameters. The plot uses the same parameters as any other biplot.
}
\value{
Returns a graphical representation of the ordinal variable on the current plot
}
\references{
Vicente-Villardon, J. L., & Sanchez, J. C. H. (2014). Logistic Biplots for Ordinal Data with an Application to Job Satisfaction of Doctorate Degree Holders in Spain. arXiv preprint arXiv:1405.0294.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
##---- Should be DIRECTLY executable !! ----
}