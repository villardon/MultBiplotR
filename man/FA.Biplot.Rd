\name{FA.Biplot}
\alias{FA.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot for Factor Analysis.
}
\description{
Biplot used as a graphical representation of Factor Analysis.
}
\usage{
FA.Biplot(X, dimension = 3, Extraction="PC", Rotation="varimax", 
         InitComunal="A1", normalize=FALSE, Scores= "Regression",  
         MethodArgs=NULL, sup.rows = NULL, sup.cols = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data Matrix
}
  \item{dimension}{
Dimension of the solution
}
  \item{Extraction}{
Method for the extraction of the factors. Can be "PC", "IPF" or "ML" ("Principal Components", "Iterated Principal Factor" or "Maximum Likelihood")
}
  \item{Rotation}{
Method for the rotation of the factors. Can be "PC", "IPF" or "ML"
}
  \item{InitComunal}{
Initial communalities for the iterated principal factor method. Can be "A1", "HSC" or "MC" ("All 1", "Highest Simple Correlation" or "Multiple Correlation")
}
  \item{normalize}{
Should the loadings be normalized
}
  \item{Scores}{
Method to calculate the Row Scores. Must be "Regression" or "Bartlett".
}
  \item{MethodArgs}{
Aditional arguments associated to the rotation method.
}
  \item{sup.rows}{
Supplementary or illustrative rows, if any.
}
  \item{sup.cols}{
Supplementary or illustrative rows, if any.
}
  \item{...}{
Additional arguments for the rotation procedure.
}
}
\details{
Biplots represent the rows and columns of a data matrix in reduced dimensions. Usually rows represent  individuals, objects or samples and columns are variables measured on them. The most classical versions can be thought as visualizations associated to Principal Components Analysis (PCA) or Factor Analysis (FA) obtained from a Singular Value Decomposition or a related method. From another point of view, Classical Biplots could be obtained from regressions and calibrations that are essentially an alternated least squares algorithm equivalent to an EM-algorithm when data are normal
This routine Calculates a biplot as a graphical representation of a classical Factor Analysis alowing for different extraction methods and different rotations.
}


\value{
An object of class "ContinuousBiplot" with the following components:

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

Gabriel, K.R.(1971): The biplot graphic display of matrices with applications to principal component analysis. Biometrika,  58, 453-467.

Gabriel, K. R. AND Zamir, S. (1979). Lower rank approximation of matrices by least squares with any choice of weights. Technometrics, 21(21):489--498, 1979.

Gabriel, K.R.(1998): Generalised Bilinear Regression. Biometrika,  85, 3, 689-700.

Gower y Hand (1996): Biplots. Chapman & Hall.

Vicente-Villardon, J. L., Galindo, M. P. and Blazquez-Zaballos, A. (2006). Logistic Biplots. Multiple Correspondence Analysis and related methods 491-509.

}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{InitialTransform}}
}
\examples{
data(Protein)
X=Protein[,3:11]
bip=FA.Biplot(X, Extraction="ML", Rotation="oblimin")
plot(bip, mode="s", margin=0.05, AddArrow=TRUE)

}

