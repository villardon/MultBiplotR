\name{plot.StatisBiplot}
\alias{plot.StatisBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a Statis Biplot Object
}
\description{
Plots a Statis Biplot Object
}
\usage{
\method{plot}{StatisBiplot}(x, A1 = 1, A2 = 2, PlotType = "Biplot", 
PlotRowTraj = FALSE, PlotVarTraj = FALSE, LabelTraj = "Begining", 
VarColorType = "ByVar", VarColors = NULL, VarLabels = NULL, 
RowColors = NULL, TableColors = NULL, RowRandomColors = FALSE, 
TypeTraj = "line", ...)
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
  \item{PlotType}{
Type of plot: Interstructure, Correlations, Contributions or Biplot
}
  \item{PlotRowTraj}{
Should the row trajectories be plotted?
}
  \item{PlotVarTraj}{
Should the variables trajectories be plotted?
}
  \item{LabelTraj}{
Where the trajecories should be labelled: Begining or End.
}
  \item{VarColorType}{
The colors for the variables should be set by table (ByTable) or by variable (ByVar)
}
  \item{VarColors}{
Colors for the variables.
}
  \item{VarLabels}{
Labels for the variables
}
  \item{RowColors}{
Colors for the rows
}
  \item{TableColors}{
Colors for each table
}
  \item{RowRandomColors}{
Use random colors for the variables.
}
  \item{TypeTraj}{
Type of trajectory to plot: Lines or stars
}
  \item{\dots}{
Aditional parameters}
}
\details{
Plots a Statis Biplot Object. The arguments of the general biplot are as in a Continuous Biplot.
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
\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
data(Chemical)
x= Chemical[,5:16]
X=Convert2ThreeWay(x,Chemical$WEEKS, columns=FALSE)
stbip=StatisBiplot(X)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
