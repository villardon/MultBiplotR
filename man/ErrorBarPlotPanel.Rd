\name{ErrorBarPlotPanel}
\alias{ErrorBarPlotPanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a panel of error bars
}
\description{
Plots a panel of error bars to compare the means of several variables in the levels of a factor using confidence intervals.
}
\usage{
ErrorBarPlotPanel(X, groups = NULL, nrows = NULL, panel = TRUE, 
GroupsTogether = TRUE, Confidence = 0.95, p.adjust.method = "None", 
UseANOVA = FALSE, Colors = "blue", Title = "Error Bar Plot", 
sort = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A matrix containing several variables
}
  \item{groups}{
A factor defining groups of individuals
}
  \item{nrows}{
Number of rows of the panel. The function calculates the number of columns needed.
}
  \item{panel}{
The plots are shown on a panel (TRUE) or in separated graphs (FALSE)
}
  \item{GroupsTogether}{
The groups appear together on the same plot
}
  \item{Confidence}{
Confidence levels for the error bars (confidence intervals)
}
  \item{p.adjust.method}{
Method for adjusting the p-value to cope with multiple comparisons.
}
  \item{UseANOVA}{
If TRUE the function uses the residual variance of the ANOVA to calculate the confidence interval. ("None", "Bonferroni" or "Sidak")
}
  \item{Colors}{
Colors to identyfy the groups
}
  \item{Title}{
Title of the graph
}
  \item{sort}{
Should short the means before plotting
}
  \item{\dots}{
Other graphical parameters
}
}

\details{
The funtion plots a panel of error bars plots to compare several groups for several variables.
}
\value{
A panel of error bars plots.
}

\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
ErrorBarPlotPanel(wine[4:9], wine$Group, UseANOVA=TRUE, Title="", sort=FALSE)
}

