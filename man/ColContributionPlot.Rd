\name{ColContributionPlot}
\alias{ColContributionPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the contributios of a biplot
}
\description{
Plots the contributios of a biplot
}
\usage{
ColContributionPlot(bip, A1 = 1, A2 = 2, Colors = NULL, Labs = NULL, 
MinQuality = 0, CorrelationScale = FALSE, ContributionScale = TRUE, 
AddSigns2Labs = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bip}{
An object of class ContinuousBiplot
}
  \item{A1}{
First dimension to plot
}
  \item{A2}{
Second dimension to plot
}
  \item{Colors}{
Colors for the variables
}
  \item{Labs}{
Labels for the variables
}
  \item{MinQuality}{
Min quality to plot
}
  \item{CorrelationScale}{
Scales for correlation
}
  \item{ContributionScale}{
Scales for contributions
}
  \item{AddSigns2Labs}{
Add the siggns of the correlations to the labels
}
  \item{\dots}{
Any other graphical parameter
}
}
\details{
Plots the contributions on a plot that sows also the sum for both axes-
}
\value{
The contribution plot
}

\author{
Jose Luis Vicente Villardon
}


\examples{
\donttest{## Simple Biplot with arrows
data(Protein)
bip=PCA.Biplot(Protein[,3:11])

# Plot of the Variable Contributions
ColContributionPlot(bip, cex=1)}
}
