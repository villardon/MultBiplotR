\name{plot.PCA.Analysis}
\alias{plot.PCA.Analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a Principal Component Analysis
}
\description{
Plots the results of a Principal Component Analysis.
}
\usage{
\method{plot}{PCA.Analysis}(x, A1 = 1, A2 = 2, CorrelationCircle = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The object with the results of a PCA
}
  \item{A1}{
Dimension for the first axis of the representation
}
  \item{A2}{
Dimension for the second axis of the representation
}
  \item{CorrelationCircle}{
Should the correlation circle be plotted? If false the scores plot is done.
}
  \item{\dots}{
Any other arguments of the function plot.ContinuousBiplot
}
}
\details{
Plots theresults of a Principal Component Analysis. The plot can be the correlation circle containing the correlations of the variables with the components or a plot of the scores of the individuals.
}
\value{
The PCA plot.
}

\author{
Jose Luis Vicente Villardon
}


\seealso{
\code{\link{plot.ContinuousBiplot}}
}
\examples{
# Not yet
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
