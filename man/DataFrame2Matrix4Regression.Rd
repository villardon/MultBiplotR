\name{DataFrame2Matrix4Regression}
\alias{DataFrame2Matrix4Regression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepares a matrix for regression from a data frame
}
\description{
Prepares a matrix for regression from a data frame
}
\usage{
DataFrame2Matrix4Regression(X, last = TRUE, Intercept = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame
}
  \item{last}{
Logical to use the last category of nominal variabless as baseline.
}
  \item{Intercept}{
Logical to tell the function if a constant must be present
}
}
\details{
Nominal variables are converted to a matrix of dummy variables for regression. 
}
\value{
A matrix ready to use as independent variables in a regression
}

\author{
Jose Luis Vicente Vilardon
}

\examples{
##---- Should be DIRECTLY executable !! ----
}

