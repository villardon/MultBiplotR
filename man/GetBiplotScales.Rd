\name{GetBiplotScales}
\alias{GetBiplotScales}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the scales for the variables on a linear biplot
}
\description{
Calculates the scales for the variables on a linear prediction biplot
There are several types of scales and values that can be shown on the graphical representation. See details.
}
\usage{
GetBiplotScales(Biplot, nticks = 4, TypeScale = "Complete", ValuesScale = "Original")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Biplot}{
Object of class PCA.Biplot
}
  \item{nticks}{
Number of ticks for the biplot axes
}
  \item{TypeScale}{
Type of scale to use : "Complete", "StdDev" or "BoxPlot"
}
  \item{ValuesScale}{
Values to show on the scale: "Original" or "Transformed"
}
}
\details{
The function calculates the points on the biplot axes where the scales should be placed.

There are three types of scales when the transformations of the raw data are made by columns:

"Complete": Covers the whole range of the variable using the number of ticks specified in "nticks". A smaller number of points could be shown if some fall outsite the range of the scatter.

"StdDev": The mean +/- 1, 2 and 3 times the standard deviation.A smaller number of points could be shown if some fall outsite the range of the scatter.

"BoxPlot": Median, 25, 75 percentiles maximum and minimum values are shown. The extremes of the interquartile range are connected with a thicker line. A smaller number of points could be shown if some fall outsite the range of the scatter.

There are two kinds of values that can be shown on the biplot axis:

"Original": The values before transformation. Only makes sense when the transformations are for each column.

"Transformed": The values after transformation, for example, after standardization. 

Although the function is public, the end used will not normally use it.
}
\value{
A list with the following components:
\item{Ticks}{A list containing the ticks for each variable}
\item{Labels}{A list containing the labels for each variable}
%% ...
}

\author{
Jose Luis Vicente Villardon
}


\examples{
data(iris)
bip=PCA.Biplot(iris[,1:4])
GetBiplotScales(bip)
}

