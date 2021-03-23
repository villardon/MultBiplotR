\name{VarBiplot}
\alias{VarBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draws a variable on a biplot
}
\description{
Draws a continuous variable on a biplot
}
\usage{
VarBiplot(bi1, bi2, b0 = 0, xmin = -3, xmax = 3, ymin = -3, ymax
                 = 3, label = "Point", mode = "a", CexPoint = 0.8,
                 PchPoint = 1, Color = "blue", ticks = c(-3, -2.5, -2,
                 -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3), ticklabels =
                 round(ticks, digits = 2), tl = 0.03, ts = "Complete",
                 Position = "Angle", AddArrow=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bi1}{
First component of the direction vector
}
  \item{bi2}{
Second component of the direction vector
}
  \item{b0}{
Constant for the regression adjusted biplots
}
  \item{xmin}{
Minimum value of the x axis
}
  \item{xmax}{
Maximum value of the x axis
}
  \item{ymin}{
Minimum value of the y axis
}
  \item{ymax}{
Maximum value of the y axis
}
  \item{label}{
Label of the variable
}
  \item{mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{CexPoint}{
Size for the symbols and labels of the variables
}
  \item{PchPoint}{
  Symbols for the variable (when represented as a point)
}
  \item{Color}{
Color for the variable
}
  \item{ticks}{
Ticks when the variable is represented as a graded scale
}
  \item{ticklabels}{
Labels for the ticks when the variable is represented as a graded scale
}
  \item{tl}{
Thick length
}
  \item{ts}{
Size of the mark in the gradedv  scale
}
  \item{Position}{
If the Position is "Angle" the label of the variable is placed using the angle of the vector
}
\item{AddArrow}{
Add an arrow to the representation of other modes of the biplot.
}
  \item{...}{
Any other graphical parameters 
}
}
\details{
See plot.PCA.Biplot
}
\value{
No value returned
}

\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
data(Protein)
bip=PCA.Biplot(Protein[,3:11])
plot(bip)
}

