\name{DimensionLabels}
\alias{DimensionLabels}

\title{
Labels for the selected dimensions in a biplot
}
\description{
Creates a character vector with labels for the dimensions of the biplot
}
\usage{
DimensionLabels(dimens, Root = "Dim")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dimens}{
Number of dimensions
}
  \item{Root}{
Root of the label
}
}
\details{
An auxiliary function to cretae labels for the dimensions. Useful to label the matrices of results
}
\value{
Returns a vector of labels
}

\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
DimensionLabels(dimens=3, Root = "Dim")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
