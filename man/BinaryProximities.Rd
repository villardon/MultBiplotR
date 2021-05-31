\name{BinaryProximities}
\alias{BinaryProximities}

\title{Proximity Measures for Binary Data}
\description{Calculation of proxymities among rows or columns of a binary data matrix or a 
data frame that will be converted into a binary data matrix. 
}
\usage{
BinaryProximities(x, y = NULL, coefficient = "Jaccard", transformation =
                 NULL, transpose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A data frame or a binary data matrix. Proximities among the rows of \code{x} will be calculated
}
  \item{y}{
Supplementary data. The proximities amond the rows of \code{x} and the rows of \code{y}  will be also calculated
}
  \item{coefficient}{
Similarity coefficient. Use the number or the name (see details)
}
  \item{transformation}{
Transformation of the similarities. Use the number or the name (see details)
}
  \item{transpose}{
Logical. If \code{TRUE}, proximities among columns are calculated
}
  \item{\dots}{
Used to provide additional parameters for the conversion of the dataframe into a binary matrix
}
}
\details{

A binary data matrix is a matrix with values 0 or 1 coding the absence or presence of several binary characters. When a data frame is provided, every variable in the data frame is converted to a binary variable using the function \code{Dataframe2BinaryMatrix}. Factors with two levels are converted directly to binary variables, factors with more than two levels are converted to a matrix with as meny columns as levels and numerical variables are converted to binary variables using a cut point that can be the median, the mean or a value provided by the user.


The following coefficients are calculated

1.- Kulezynski =  a/(b + c)

2.- Russell_and_Rao =  a/(a + b + c+d)

3.- Jaccard = a/(a + b + c)

4.- Simple_Matching = (a + d)/(a + b + c + d)

5.- Anderberg = a/(a + 2 * (b + c))

6.- Rogers_and_Tanimoto = (a + d)/(a + 2 * (b + c) + d)

7.- Sorensen_Dice_and_Czekanowski = a/(a + 0.5 * (b + c))

8.- Sneath_and_Sokal = (a + d)/(a + 0.5 * (b + c) + d)

9.- Hamman = (a - (b + c) + d)/(a + b + c + d)

10.- Kulezynski = 0.5 * ((a/(a + b)) + (a/(a + c)))

11.- Anderberg2 = 0.25 * (a/(a + b) + a/(a + c) + d/(c + d) + d/(b + d))

12.- Ochiai = a/sqrt((a + b) * (a + c))

13.- S13 = (a * d)/sqrt((a + b) * (a + c) * (d + b) * (d + c))

14.- Pearson_phi = (a * d - b * c)/sqrt((a + b) * (a + c) * (d + b) * (d + c))

15.- Yule = (a * d - b * c)/(a * d + b * c)

The following transformations of the similarity3 are calculated

1.- `Identity`   dis=sim

2.- `1-S`   dis=1-sim

3.- `sqrt(1-S)`  dis = sqrt(1 - sim)

4.- `-log(s)`  dis=-1*log(sim)

5.- `1/S-1` dis=1/sim -1 

6.- `sqrt(2(1-S))`  dis== sqrt(2*(1 - sim))

7.- `1-(S+1)/2` dis=1-(sim+1)/2

8.- `1-abs(S)`  dis=1-abs(sim)

9.- `1/(S+1)`  dis=1/(sim)+1

Note that, after transformation the similarities are converted to distances except for "Identity".
Not all the transformations are suitable for all the coefficients. Use them at your own risk.
The default values are admissible combinations.
  
}
\value{
An object of class \code{proximities}.This has components:
\item{TypeData }{Binary, Continuous or Mixed. Binary in this case.}
\item{Coefficient }{Coefficient used to calculate the proximities}
\item{Transformation }{Transformation used to calculate the proximities}
\item{Data }{Data used to calculate the proximities}
\item{SupData }{Supplementary Data, if any}
\item{Proximities }{Proximities among rows of \code{x}. May be similarities or dissimilarities depending on the transformation}
\item{SupProximities }{Proximities among rows of \code{x} and \code{y}.}
}

\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley
}
\author{
Jose Luis Vicente-Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BinaryDistances}}, \code{\link{Dataframe2BinaryMatrix}}
}
\examples{
\donttest{data(spiders)
D=BinaryProximities(spiders, coefficient="Jaccard", transformation="sqrt(1-S)")
D2=BinaryProximities(spiders, coefficient=3, transformation=3)}
}

