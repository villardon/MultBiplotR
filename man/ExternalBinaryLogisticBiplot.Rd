\name{ExternalBinaryLogisticBiplot}
\alias{ExternalBinaryLogisticBiplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
External Logistic Biplot for binary Data
}
\description{
Fits an External Logistic Biplot to the results of a Principal Coordinates Analysis obtained from binary data.

}
\usage{
ExternalBinaryLogisticBiplot(Pco, IncludeConst=TRUE,  penalization=0.2, freq=NULL, 
tolerance = 1e-05, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Pco}{
An object of class "Principal.Coordinates"
}
  \item{IncludeConst}{
Should the logistic fit include the constant term?
}
  \item{penalization}{
Penalization for the ridge regression   
}

  \item{freq}{
frequencies for each observation or pattern (usually 1)
}

  \item{tolerance}{
Tolerance for convergence
}

\item{maxiter}{
Maximum number of iterations
}

}
\details{
Let \eqn{{\bf{X}}}  be the matrix of binary data scored as present or absent (1 or 0), in which the rows correspond to n individuals or entries (for example, genotypes) and the columns to p binary characters (for example alleles or bands), let  \eqn{{\bf{S}} = ({s_{ij}})} be a matrix containing the similarities among rows, obtained from the binary data matrix , and let \eqn{\Delta  = ({\delta _{ij}})}  be the corresponding dissimilarity/distance matrix, taking for example \eqn{{\delta _{ij}} = \sqrt {1 - {s_{ij}}}}.
Despite the fact that, in Cluster Analysis and Principal Coordinates Analysis, interpretation of the variables responsible for grouping or ordination is not straightforward, those methods are normally used to classify individual in which binary variables have been measured. 
we use a combination of Principal Coordinates Analysis (PCoA), Cluster Analysis (CA) and External Logistic Regression (ELB), as a better way to interpret the binary variables associated to the classification of genotypes. The combination of three standard techniques with some new ideas about the geometry of the procedures, allows to construct a External Logistic Regression (ELB), that helps the interpretation of the variables responsible for the classification or ordination.
Suppose we have obtained an euclidean configuration \eqn{{\bf{Y}}} obtained from the Principal Coordinates (PCoA) of the similarity matrix.
To search for the variables associated to the ordination obtained in PCoA, we can look for the directions in the ordination diagram that better predict the probability of presence of each allele.
More formally, if we defined \eqn{{\pi _{ij}} = E({x_{ij}})= {\textstyle{1 \over {1 + \exp ( - ({b_{j0}} + \sum\limits_{s = 1}^k {{b_{js}}{y_{is}}} ))}}}} as the expected probability that the allele j be present at genotype for a genotype with coordinates \eqn{y_{is}}  (i=1, ...,n; s=1, ..., k) on the ordination diagram, as
where bjs ( j=1,..., p) are the logistic regression coefficients that correspond to the jth variable (alleles or bands) in the sth dimension. The model is a generalized linear model having the logit as a link function.
where    and , y's and b's define a biplot in logit scale. This is called External Logistic Biplot because the coordinates of the genotypes are calculated in an external procedure (PCoA). Given that the y's are known from PCoA, obtaining the b´s is equivalent to performing a logistic regression using the j-th column of X as a response variable and the columns of y as regressors.
}
\value{
An object of class \code{External.Binary.Logistic.Biplot} with the fields of the \code{Principal.Coordinates} object with the following fields added.
\item{ColumnParameters }{Parameters resulting from fitting a logistic regression to each column of the original binary data matrix}
\item{VarInfo}{Information of the fit for each variable}
\item{VarInfo$Deviances}{A vector with the deviances of each variable calculated as the difference with the null model}
\item{VarInfo$Dfs}{A vector with degrees of freedom for each variable}
\item{VarInfo$pvalues}{A vector with the p values each variable}
\item{VarInfo$Nagelkerke }{A vector with the Nagelkerke pseudo R-squared for each variable}
\item{VarInfo$PercentsCorrec}{A vector with the percentage of correct classifications for each variable}
\item{DevianceTotal}{Total Deviance as the difference with the null model}
\item{p}{p value for the complete representation}
\item{TotalPercent}{Total percentage of correct classification}

}
\references{
Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.

Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Análisis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(spiders)
x2=Dataframe2BinaryMatrix(spiders)
colnames(x2)=colnames(spiders)
dist=BinaryProximities(x2)
pco=PrincipalCoordinates(dist)
pcobip=ExternalBinaryLogisticBiplot(pco)
}

