\name{BasicDescription}
\alias{BasicDescription}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Basic descriptive sataistics
}
\description{
Basic descriptive sataistics of several variables by the categories of a factor.
}
\usage{
BasicDescription(X, groups = NULL, SortByGroups = FALSE, na.rm = FALSE, Intervals = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame or a matrix containing several numerical variables
}
  \item{groups}{
A factor with the groupings
}
  \item{SortByGroups}{
Sorting by groups
}
  \item{na.rm}{
a logical value indicating whether NA values should be stripped before the computation proceeds.
}
  \item{Intervals}{
Should the confidence intervals be calculated?
}
}
\details{
Basic descriptive sataistics of several variables by the categories of a factor.
}
\value{
A list with the description of each variable.
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
BasicDescription(wine[,4:8], groups = wine$Origin)
}
