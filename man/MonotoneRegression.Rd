\name{MonotoneRegression}
\alias{MonotoneRegression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Isotonic Regression (Weighted Monotone Regression)
}
\description{
Performs weighted isotonic (monotone) regression using the non-negative weights in w. The function is a direct translation of the matlab function lsqisotonic.
}
\usage{
MonotoneRegression(x, y, w = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The independent variable vector
}
  \item{y}{
The dependent variable vector
}
  \item{w}{
A vector of weigths
}
}
\details{
YHAT = MonotoneRegression(X,Y) returns a vector of values that minimize the
sum of squares (Y - YHAT).^2 under the monotonicity constraint that
X(I) > X(J) => YHAT(I) >= YHAT(J), i.e., the values in YHAT are
monotonically non-decreasing with respect to X (sometimes referred
to as "weak monotonicity").  LSQISOTONIC uses the "pool adjacent
violators" algorithm.

If X(I) == X(J), then YHAT(I) may be <, ==, or > YHAT(J) (sometimes
referred to as the "primary approach").  If ties do occur in X, a plot
of YHAT vs. X may appear to be non-monotonic at those points.  In fact,
the above monotonicity constraint is not violated, and a reordering
within each group of ties, by ascending YHAT, will produce the desired
appearance in the plot.
}
\value{
The fitted values after the monotone regression
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kruskal, J.B. (1964) "Nonmetric multidimensional scaling: a numerical method", Psychometrika 29:115-129.

Cox, R.F. and Cox, M.A.A. (1994) Multidimensional Scaling, Chapman&Hall.
}
\author{
Jose L. Vicente Villardon (from a matlab functiom)
}
\note{
The function is a direct translation of the matlab function lsqisotonic.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Used inside MDS

}

