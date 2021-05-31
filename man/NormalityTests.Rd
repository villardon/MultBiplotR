\name{NormalityTests}
\alias{NormalityTests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normality tests
}
\description{
Normality tests foor the columns of a matrix  and a grouping variable.
}
\usage{
NormalityTests(X, groups = NULL, plot = FALSE, SortByGroups = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame or a matrix containing several numerical variables
}
  \item{groups}{
A factor with the groups
}
  \item{plot}{
If TRUE the qqnorm plots are shown
}
  \item{SortByGroups}{
Should the results be sorted by groups?
}
}
\details{
Normality tests foor the columns of a matrix  and a grouping variable.
}
\value{
The normality tests and the plots
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
NormalityTests(wine[,4:6], groups = wine$Origin, plot=TRUE)
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
