\name{NIPALS.Biplot}
\alias{NIPALS.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot using the NIPALS algorithm
}
\description{
Biplot using the NIPALS algorithm including a truncated and a sparse version.
}
\usage{
NIPALS.Biplot(X, alpha = 1, dimension = 3, Scaling = 5, 
Type = "Regular", grouping = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The data matrix
}
  \item{alpha}{
A number between 0 and 1. 0 for GH-Biplot, 1 for JK-Biplot and 0.5 for SQRT-Biplot. Use 2 or any other value not in the interval [0,1] for HJ-Biplot.
}
  \item{dimension}{
Dimension of the solution
}
  \item{Scaling}{
Transformation of the original data. See InitialTransform for available transformations.
}
  \item{Type}{
Type of biplot (Regular, Truncated or Sparse)
}
  \item{grouping}{
Grouping fartor when the scaling is made with the within groups variability
}
  \item{\dots}{
Aditional arguments for the different types of biplots.
}
}
\details{
Biplot using the NIPALS algorithm including a truncated and a sparse version.
}
\value{
An object of class ContinuousBiplot with the following components:


\item{Title }{A general title}
\item{Type }{NIPALS}
\item{call }{call}
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
Wold, H. (1966). Estimation of principal components and related models by iterative least squares. Multivariate analysis. ACEDEMIC PRESS. 391-420.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
bip1=NIPALS.Biplot(wine[,4:21], Type="Sparse", lambda=0.15)
plot(bip1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
