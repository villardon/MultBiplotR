\name{TransformIni}
\alias{TransformIni}
 
\title{
Initial transformation of a data matrix
}
\description{
Initial transformation of data before the construction of a biplot. (or any other technique)
}
\usage{
TransformIni(X, InitTransform = "None", transform = "Standardize columns")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Original Raw Data Matrix
}
  \item{InitTransform}{
Initial transform of the data (usually logarithm)
}
  \item{transform}{
Transformation to use. See details.
}
}
\details{
Possible Transformations are:

1.- "Raw Data": When no transformation is required.

2.- "Substract the global mean": Eliminate an eefect common to all the observations

3.- "Double centering" : Interaction residuals. When all the elements of the table are comparable. Useful for AMMI models.

4.- "Column centering": Remove the column means.

5.- "Standardize columns":  Remove the column means and divide by its standard deviation.

6.- "Row centering": Remove the row means.

7.- "Standardize rows": Divide each row by its standard deviation.

8.- "Divide by the column means and center": The resulting dispersion is the coefficient of variation.

9.- "Normalized residuals from independence" for a contingency table.

The transformation can be provided to the function by using the string beetwen the quotes or just the associated number.

The supplementary rows and columns are not used to calculate the parameters (means, standard deviations, etc). Some of the transformations are not compatible with supplementary data. 
}
\value{

\item{X }{Transformed data matrix}
 
}
\references{
M. J. Baxter (1995) Standardization and Transformation in Principal Component Analysis, with Applications to Archaeometry. 
Journal of the Royal Statistical Society. Series C (Applied Statistics). Vol. 44, No. 4 (1995) , pp. 513-527

Kroonenberg, P. M. (1983). Three-mode principal component analysis: Theory and applications (Vol. 2). DSWO press. (Chapter 6)
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(iris)
x=as.matrix(iris[,1:4])
x=TransformIni(x, transform=4)
x
}

