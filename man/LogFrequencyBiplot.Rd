\name{LogFrequencyBiplot}
\alias{LogFrequencyBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Biplot for a table of frequencies 
}
\description{
Biplot for the logarithms of the frequencies of a contingency table using the frequencies as weights.
}
\usage{
LogFrequencyBiplot(x, Scaling = 2, logoffset = 1, freqoffset = logoffset, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The frequency table to be biplotted
}
  \item{Scaling}{
Transformation of the matrix after the logarithms
}
  \item{logoffset}{
Constant to add to the frequencies before calculating the logarithms. This is to avoid calculating the logaritm of zero, so, a covenient value for this argument is 1.
}
  \item{freqoffset}{
Constant to add to the frequencies before calculating the weigths. This is usually the same as the offset used to add to the frequencies but may be different when we do not want the frequencies zero to influence the biplot, i. e., we want zero weigths.
}
  \item{\dots}{
Any other parameter for the CrissCross procedure.
}
}
\details{
Biplot for the logarithms of the frequencies of a contingency table using the frequencies as weigths.
}
\value{
An object of class .Biplot" with the following components:


\item{Title }{A general title}
\item{Non_Scaled_Data }{Original Data Matrix}
\item{Means }{Means of the original Variables}
\item{Medians }{Medians of the original Variables}
\item{Deviations }{Standard Deviations of the original Variables}
\item{Minima }{Minima of the original Variables}
\item{Maxima }{Maxima of the original Variables}
\item{P25 }{25 Percentile of the original Variables}
\item{P75 }{75 Percentile of the original Variables}
\item{Gmean }{Global mean of the complete matrix}
\item{Sup.Rows }{Supplementary rows (Non Transformed)}
\item{Sup.Cols }{Supplementary columns (Non Transformed)}
\item{Scaled_Data }{Transformed Data}
\item{Scaled_Sup.Rows}{Supplementary rows (Transformed)}
\item{Scaled_Sup.Cols}{Supplementary columns (Transformed)}
\item{n}{Number of Rows}
\item{p}{Number of Columns}
\item{nrowsSup}{Number of Supplementary Rows}
\item{ncolsSup}{Number of Supplementary Columns}
\item{dim}{Dimension of the Biplot}
\item{EigenValues}{Eigenvalues}
\item{Inertia}{Explained variance (Inertia)}
\item{CumInertia}{Cumulative Explained variance (Inertia)}
\item{EV}{EigenVectors}
\item{Structure}{Correlations of the Principal Components and the Variables}
\item{RowCoordinates}{Coordinates for the rows, including the supplementary}
\item{ColCoordinates}{Coordinates for the columns, including the supplementary}
\item{RowContributions}{Contributions for the rows, including the supplementary}
\item{ColContributions}{Contributions for the columns, including the supplementary}
\item{Scale_Factor}{Scale factor for the traditional plot with points and arrows. The row coordinates are multiplied and the column coordinates divided by that scale factor. The look of the plot is better without changing the inner product. For the HJ-Biplot the scale factor is 1.}
}
\references{
Gabriel, K. R., Galindo, M. P. & Vicente-Villardon, J. L. (1995) Use of Biplots to Diagnose Independence Models in Three-Way Contingency Tables. in: M. Greenacre & J. Blasius. eds. Visualization of Categorical Data. Academis Press. London.

GABRIEL, K.R. and ZAMIR, S. (1979). Lower rank approximation of matrices by least squares with any choice of weights. Technometrics, 21: 489-498.
}
\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{CrissCross}}, ~~~
}
\examples{
data(smoking)
logbip=LogFrequencyBiplot(smoking, Scaling=1, logoffset=0, freqoffset=0)
}

