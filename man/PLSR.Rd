\name{PLSR}
\alias{PLSR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Partial Least Squares Regression
}
\description{
Partial Least Squares Regression for numerical variables.
}
\usage{
PLSR(Y, X, S = 2, InitTransform = 5, grouping = NULL, 
centerY = TRUE, scaleY = TRUE, tolerance = 5e-06, 
maxiter = 100, show = FALSE, Validation = NULL, nB = 500)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of Dependent Variables
}
  \item{X}{
Matrix of Independent Variables
}
  \item{S}{
Dimension of the solution
}
  \item{InitTransform}{
Initial transformation of the independent variables.
}
  \item{grouping}{
Fator when the init transformation is the standardization with the within groups deviation.
}
  \item{centerY}{
Should the dependent variables be centered?
}
  \item{scaleY}{
Should the dependent variables be standadized?
}
  \item{tolerance}{
Tolerance for the algorithm
}
  \item{maxiter}{
Maximum number of iterations
}
  \item{show}{
Show the progress of the algorithm?
}
  \item{Validation}{
Validation (None, Cross, Bootstrap)
}
  \item{nB}{
number of samples for the bottstrap validation
}
}
\details{
Partial Least Squares Regression for numerical variables.
}
\value{
An object of class plsr with fiends

\item{Method}{PLSR}
\item{X}{The X matrix}
\item{Y}{The Y matrix}
\item{centerY}{Is the Y matrix centered}
\item{scaleY}{Is the Y matrix scaled}
\item{Initial_Transformation}{Initial transformation of the Y matrix}
\item{ScaledX}{Transformed X matrix}
\item{ScaledY}{Transformed Y matrix}
\item{Intercept}{Intercept of the model}
\item{XScores}{Scores for the individals from the X matrix}
\item{XWeights}{Weigths for the X set}
\item{XLoadings}{Loadings for the X set}
\item{YScores}{Scores for the individals from the Y matrix}
\item{YWeights}{Weigths for the Y set}
\item{YLoadings}{Loadings for the Y set}
\item{RegParameters}{Final Regression Parameters}
\item{ExpectedY}{Expected values of Y}
\item{R2}{R-squared}
\item{XStructure}{Relation of the X variables with its structure}
\item{YStructure}{Relation of the Y variables with its structure}
\item{YXStructure}{Relation of the Y variables with the X components}
}
\references{
H. Abdi, Partial least squares regression and projection on latent structure regression (PLS regression), WIREs Comput. Stat. 2 (2010), pp. 97-106.
}
\author{
Jose Luis Vicente Villardon
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Biplot.PLSR}}
}
\examples{
X=as.matrix(wine[,4:21])
y=as.numeric(wine[,2])-1
mifit=PLSR(y,X, Validation="None")
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
