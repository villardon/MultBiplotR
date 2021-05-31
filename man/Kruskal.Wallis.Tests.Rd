\name{Kruskal.Wallis.Tests}
\alias{Kruskal.Wallis.Tests}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kruskal Wallis Tests
}
\description{
Kruskal Wallis Tests for a matrix of continuous variables and a grouping factor.
}
\usage{
Kruskal.Wallis.Tests(X, groups, posthoc = "none", alternative = "two.sided", digits = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The matrix of continuous variables
}
  \item{groups}{
The factor with the groups
}
  \item{posthoc}{
Method used for multipe comparisons in the Dunn test
}
  \item{alternative}{
Kind of alternative hypothesis
}
  \item{digits}{
number of digitd for he output
}
}
\details{
Kruskal Wallis Tests for a matrix of continuous variables and a grouping factor, including the Dunn test for multiple comparisons.
}
\value{
the organized output.
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
Kruskal.Wallis.Tests(wine[,4:7], wine$Group, posthoc = "bonferroni")
}

