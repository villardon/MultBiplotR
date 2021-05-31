\name{Canonical.Variate.Analysis}
\alias{Canonical.Variate.Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Biplot representation of a Canonical Variate Analysis or a Manova (Canonical-Biplot or MANOVA-Biplot)
}
\description{
Calculates a canonical biplot with confidence regions for the means.
}
\usage{
Canonical.Variate.Analysis(X, group, InitialTransform = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data matrix
}
  \item{group}{
A factor containing the groups
}
\item{InitialTransform}{
Initial transformation of the data matrix
}
}
\details{
The Biplot method (Gabriel, 1971; Galindo, 1986; Gower and Hand, 1996) is becoming one of the most popular techniques for analysing multivariate data. Biplot methods are techniques for simultaneous representation of the \eqn{n} rows and \eqn{n} columns of a data matrix \eqn{\bf{X}}, in reduced dimensions, where the rows represent individuals, objects or samples and the columns the variables measured on them. Classical Biplot methods are a graphical representation of a Principal Components Analysis (PCA) that it is used to obtain linear combinations that successively maximize the total variability. 
PCA is not considered an appropriate approach where there is known a priori group structure in the data. The most general methodology for discrimination among groups, using multiple observed variables, is Canonical Variate Analysis (CVA). CVA allows us to derive linear combinations that successively maximize the ratio of "between-groups"" to "pooled within-group" sample variance. Several authors propose a Biplot representation for CVA called Canonical Biplot (CB) (Vicente-Villardon, 1992 and Gower & Hand, 1996) when it is oriented to the discrimination between groups or MANOVA-Biplot Gabriel (1972, 1995) when the aim is to study the variables responsible for the discrimination. The main advantage of the Biplot version of the technique is that it is possible not only to establish the differences between groups but also to characterise the variables responsible for them. The methodology is not yet widely used mainly because it is still not available in the major statistical packages.
Amaro, Vicente-Villardon & Galindo (2004) extend the methodology for two-way designs and propose confidence circles based on univariate  and multivariate tests to perform post-hoc analysis of each variable.
}
\value{
An object of class "Canonical.Biplot"
}
\references{
Amaro, I. R., Vicente-Villardon, J. L., & Galindo-Villardon, M. P. (2004). Manova Biplot para arreglos de tratamientos con dos factores basado en modelos lineales generales multivariantes. Interciencia, 29(1), 26-32.

Vicente-Villardón, J. L. (1992). Una alternativa a las técnicas factoriales clásicas basada en una generalización de los métodos Biplot (Doctoral dissertation, Tesis. Universidad de Salamanca. España. 248 pp.[Links]).

Gabriel KR (1971) The biplot graphic display of matrices with application to principal component analysis. Biometrika 58(3):453-467.

Gabriel, K. R. (1995). MANOVA biplots for two-way contingency tables. WJ Krzanowski (Ed.), Recent advances in descriptive multivariate analysis, Oxford University Press, Toronto. 227-268.

Galindo Villardon, M. (1986). Una alternativa de representacion simultanea: HJ-Biplot. Qüestiió. 1986, vol. 10, núm. 1.

Gower y Hand (1996): Biplots. Chapman & Hall.

Varas, M. J., Vicente-Tavera, S., Molina, E., & Vicente-Villardon, J. L. (2005). Role of canonical biplot method in the study of building stones: an example from Spanish monumental heritage. Environmetrics, 16(4), 405-419.

Santana, M. A., Romay, G., Matehus, J., Villardon, J. L., & Demey, J. R. (2009). simple and low-cost strategy for micropropagation of cassava (Manihot esculenta Crantz). African Journal of Biotechnology, 8(16).

}
\author{
Jose Luis Vicente Villardon
}

\examples{
\donttest{data(wine)
X=wine[,4:21]
canbip=CanonicalBiplot(X, group=wine$Group)
plot(canbip, mode="s")}
}
