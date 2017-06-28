\name{MultBiplotR-package}
\alias{MultBiplotR-package}
\alias{MultBiplot}
\docType{package}
\title{
Multivariate Analysis using Biplots 
}
\description{
Classical PCA biplot with aditional features as non-standard data transformations, scales for the variables, together with many graphical aids as sizes or colors of the points according to their qualities of representation or predictiveness. The package includes also Alternating Least Squares (ALS) or Criss-Cross procedures for the calculation of the reduced rank approximation that can deal with missing data, differencial weights for each element of the data matrix or even ronust versions of the procedure. 

This is part of a bigger project called MULTBIPLOT that contains many other biplot techniques and is a translation to R of the package MULBIPLOT programmed in MATLAB. A GUI for the package is also in preparation.
}
\details{
\tabular{ll}{
Package: \tab MultBiplot\cr
Type: \tab Package\cr
Version: \tab 0.1.00\cr
Date: \tab 2015-01-14\cr
License: \tab GPL(>=2)\cr
}

}

\author{
Jose Luis Vicente Villardon
Maintainer: Jose Luis Vicente Villardon <villardon@usal.es>
}
\references{

Vicente-Villardon, J.L. (2010). MULTBIPLOT: A package for Multivariate Analysis using Biplots. Departamento de Estadistica. Universidad de  Salamanca. (http://biplot.usal.es/ClassicalBiplot/index.html).

Vicente-Villardon, J. L. (1992). Una alternativa a las técnicas factoriales clasicas basada en una generalización de los metodos Biplot (Doctoral dissertation, Tesis. Universidad de Salamanca. España. 248 pp.[Links]).

Gabriel KR (1971) The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3):453-467

Gabriel KR (1998) Generalised bilinear regresion, J. L. (1998). Use of biplots to diagnose independence models in three-way contingency tables. Visualization of Categorical Data. Academic Press. London.

Gabriel, K. R. (2002). Le biplot-outil d'exploration de donnes multidimensionnelles. Journal de la Societe francaise de statistique, 143(3-4).

Gabriel KR, Zamir S (1979) Lower rank approximation of matrices by least squares with any choice of weights. Technometrics 21(4):489-498.

Gower J, Hand D (1996) Biplots. Monographs on statistics and applied probability. 54. London: Chapman and Hall., 277 pp.

Galindo Villardon, M. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Qüestiió. 1986, vol. 10, núm. 1.

Demey J, Vicente-Villardon JL, Galindo MP, Zambrano A (2008) Identifying molecular markers associated with classification of genotypes using external logistic biplots. Bioinformatics 24(24):2832-2838.

Vicente-Villardon JL, Galindo MP, Blazquez-Zaballos A (2006) Logistic biplots. Multiple Correspondence Analysis and related methods pp 491-509.

Santos, C., Munoz, S. S., Gutierrez, Y., Hebrero, E., Vicente, J. L., Galindo, P., \& Rivas, J. C. (1991). Characterization of young red wines by application of HJ biplot analysis to anthocyanin profiles. Journal of Agricultural and food chemistry, 39(6), 1086-1090.

Rivas-Gonzalo, J. C., Gutierrez, Y., Polanco, A. M., Hebrero, E., Vicente, J. L., Galindo, P., \& Santos-Buelga, C. (1993). Biplot analysis applied to enological parameters in the geographical classification of young red wines. American journal of enology and viticulture, 44(3), 302-308.
}

\keyword{ package }

\examples{
data(iris)
bip=PCA.Biplot(iris[,1:4])
plot(bip)
}

