\name{RidgeMultinomialLogisticFit}
\alias{RidgeMultinomialLogisticFit}
\title{
Multinomial logistic regression with ridge penalization
}
\description{
This function does a logistic regression between a dependent variable y and 
  some independent variables x, and solves the separation problem in this type
  of regression using ridge regression and penalization.
}
\usage{
RidgeMultinomialLogisticFit(y, x, penalization = 0.2, 
tol = 1e-04, maxiter = 200, show = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Dependent variable.
}
  \item{x}{
A matrix with the independent variables.
}
  \item{penalization}{
Penalization used in the diagonal matrix to avoid singularities.
}
  \item{tol}{
Tolerance for the iterations.
}
  \item{maxiter}{
Maximum number of iterations.
}
  \item{show}{
Should the iteration history be printed?.
}
}
\details{
The problem of the existence of the estimators in logistic regression can be seen in Albert (1984), a solution for the binary case, based on the Firth's method, Firth (1993) is proposed by Heinze(2002). The extension to nominal logistic model was made by Bull (2002). All the procedures were initially developed to remove the bias but work well to avoid the problem of separation. Here we have chosen a simpler solution based on ridge estimators for logistic regression Cessie(1992).  

Rather than maximizing \eqn{{L_j}(\left. {\bf{G}} \right|{{\bf{b}}_{j0}},{{\bf{B}}_j})}{L_j (G | b_j0 , B_j)} we maximize

\deqn{{{L_j}(\left. {\bf{G}} \right|{{\bf{b}}_{j0}},{{\bf{B}}_j})} - \lambda \left( {\left\| {{{\bf{b}}_{j0}}} \right\| + \left\| {{{\bf{B}}_j}} \right\|} \right)}{L_j(G|b_jo,B_j)-\lambda(||b_j0|| + ||B_j||)}
 
Changing the values of \eqn{\lambda} we obtain slightly different solutions not affected by the separation problem.
}
\value{
An object of class \code{"rmlr"} with components

 \item{fitted }{Matrix with the fitted probabilities}
 \item{cov }{Covariance matrix among the estimates}
 \item{Y }{Indicator matrix for the dependent variable}
 \item{beta }{Estimated coefficients for the multinomial logistic regression}
 \item{stderr }{Standard error of the estimates}
 \item{logLik }{Logarithm of the likelihood}
 \item{Deviance }{Deviance of the model}
 \item{AIC }{Akaike information criterion indicator} 
 \item{BIC }{Bayesian information criterion indicator} 
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
\examples{

  # No examples yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ridge}
\keyword{logistic}
\keyword{models}
