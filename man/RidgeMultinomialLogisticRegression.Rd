\name{RidgeMultinomialLogisticRegression}
\alias{RidgeMultinomialLogisticRegression}
\title{
Ridge Multinomial Logistic Regression
}
\description{
Function that calculates an object with the fitted multinomial logistic regression for a nominal variable. It compares with the null model, so that we will be able to compare which model fits better the variable.
}
\usage{
RidgeMultinomialLogisticRegression(formula, data, penalization = 0.2,
cte = TRUE, tol = 1e-04, maxiter = 200, showIter = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
The usual formula notation (or the dependent variable)
}
  \item{data}{
The dataframe used by the formula. (or a matrix with the independent variables).
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
  \item{cte}{
Should the model have a constant? 
}
  \item{tol}{
Value to stop the process of iterations.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{showIter}{
Should the iteration history be printed?.
}
}

\value{
An object that has the following components:

 \item{fitted }{Matrix with the fitted probabilities}
 \item{cov }{Covariance matrix among the estimates}
 \item{Y }{Indicator matrix for the dependent variable}
 \item{beta }{Estimated coefficients for the multinomial logistic regression}
 \item{stderr }{Standard error of the estimates}
 \item{logLik }{Logarithm of the likelihood}
 \item{Deviance }{Deviance of the model}
 \item{AIC }{Akaike information criterion indicator} 
 \item{BIC }{Bayesian information criterion indicator} 
 \item{NullDeviance }{Deviance of the null model}
 \item{Difference }{Difference between the two deviance values}
 \item{df }{Degrees of freedom}
 \item{p }{p-value asociated to the chi-squared estimate}
 \item{CoxSnell }{Cox and Snell pseudo R squared}
 \item{Nagelkerke }{Nagelkerke pseudo R squared}
 \item{MacFaden }{MacFaden pseudo R squared}
  \item{Table }{Cross classification of observed and predicted responses}
 \item{PercentCorrect }{Percentage of correct classifications}

}
\references{
  Albert,A. & Anderson,J.A. (1984),\emph{On the existence of maximum likelihood estimates in logistic regression models}, Biometrika 71(1), 1--10.
  
  Bull, S.B., Mak, C. & Greenwood, C.M. (2002), \emph{A modified score function for multinomial logistic regression}, Computational Statistics and dada Analysis 39, 57--74.
  
  Firth, D.(1993), \emph{Bias reduction of maximum likelihood estimates}, Biometrika 80(1), 27--38
  
  Heinze, G. & Schemper, M. (2002), \emph{A solution to the problem of separation in logistic regression}, Statistics in Medicine 21, 2109--2419
  
  Le Cessie, S. & Van Houwelingen, J. (1992), \emph{Ridge estimators in logistic regression}, Applied Statistics 41(1), 191--201.
}
\author{
Jose Luis Vicente-Villardon
}
\seealso{
\code{\link{RidgeMultinomialLogisticFit}}
}
\examples{
  
  data(Protein)
  y=Protein[[2]]
  X=Protein[,c(3,11)]
  rmlr = RidgeMultinomialLogisticRegression(y,X,penalization=0.0)
  summary(rmlr)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{logistic}% __ONLY ONE__ keyword per line
