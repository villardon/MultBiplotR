\name{Levene.Tests}
\alias{Levene.Tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Levene Tests 
}
\description{
Levene Tests  for a matrix of continuous variables and a grouping factor.
}
\usage{
Levene.Tests(X, groups = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The matrix of continuous variables
}
  \item{groups}{
The factor with the groups
}
}
\details{
Levene Tests  for a matrix of continuous variables and a grouping factor.
}
\value{
The organized output
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
Levene.Tests(wine[,4:7], wine$Group)
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
