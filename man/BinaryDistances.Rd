\name{BinaryDistances}
\alias{BinaryDistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binary Distances }
\description{
Calculates distances among rows of a binary data matrix or among the rows of two binary matrices. 
The end user will use BinaryProximities rather than this function.
Input must be a matrix with 0 or 1 values.
}
\usage{
BinaryDistances(x, y = NULL, coefficient= "Simple_Matching", transformation="sqrt(1-S)")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Main binary data matrix. Distances among rows are calculated if y=NULL.
}
  \item{y}{
Second binary data matrix. If not NULL the distances among the rows of x and y are calculated 
}
  \item{coefficient}{
Similarity coefficient. Use the name (see details)
}
  \item{transformation}{
Transformation of the similarities. Use the name (see details)
}
}
\details{
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

  
}
\value{
An object of class \code{proximities}.This has components:
\item{comp1 }{Description of 'comp1'}

}
\references{
Gower, J. C. (2006) Similarity dissimilarity and Distance, measures of. Encyclopedia of Statistical Sciences. 2nd. ed. Volume 12. Wiley
}
\author{
Jose Luis Vicente-Villardon
}

\seealso{
\code{\link{PrincipalCoordinates}}
}
\examples{
data(spiders)

}

