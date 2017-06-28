\name{OrdVarCoordinates}
\alias{OrdVarCoordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coordinates of an ordinal variable on the biplot.
}
\description{
Coordinates of an ordinal variable on the biplot.
}
\usage{
OrdVarCoordinates(tr, b = c(1, 1), inf = -12, sup = 12, step = 0.01,
                 plotresponse = FALSE, label = "Item", labx = "z", laby
                 = "Probability", catnames = NULL, Legend = TRUE,
                 LegendPos = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tr}{
A vector containing the thresholds of the model, that is, the constatn for each category of the ordinal variable
}
  \item{b}{
Vector containing the common slopes for all categories of the ordinal variable
}
  \item{inf}{
The inferior limit of the values to be sampled on the biplot axis (it depends on the scale of the biplot).
}
  \item{sup}{
The superior limit of the values to be sampled on the biplot axis (it depends on the scale of the biplot).
}
  \item{step}{
Increment (step) of the squence
}
  \item{plotresponse}{
Should the item be plotted
}
  \item{label}{
Label of the item.
}
  \item{labx}{
Label for the X axis in the summary of the item.
}
  \item{laby}{
Label for the Y axis in the summary of the item.
}
  \item{catnames}{
Names of the categories.
}
  \item{Legend}{
Should a legend be plotted
}
  \item{LegendPos}{
Position of the legend.
}
}
\details{
The function calculates the coordinates of the points that define the separation among the categories of an ordinal variable projected onto an ordinal logistic biplot. 
}
\value{
An object of class OrdVarCoord
\item{z }{Values of the cut points on the scale of the biplot axis (not used)}
\item{points }{The points for the marks to be represented on the biplot.}
\item{labels }{The labels for the points}
\item{hidden }{Are there any hidden categories? (Categories whose probability is never hier than the probabilities of the rest)}
\item{cathidden}{Number of the hidden cateories}

}
\references{
Vicente-Villardon, J. L., & Sanchez, J. C. H. (2014). Logistic Biplots for Ordinal Data with an Application to Job Satisfaction of Doctorate Degree Holders in Spain. arXiv preprint arXiv:1405.0294.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
# No examples
}

