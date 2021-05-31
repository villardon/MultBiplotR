\name{HistogramPanel}
\alias{HistogramPanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Panel of histograms
}
\description{
Panel of histograms for a set of numerical variables.
}
\usage{
HistogramPanel(X, nrows = NULL, separated = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The matrix of continuous variables
}
  \item{nrows}{
Number of rows of the panel.
}
  \item{separated}{
Should the plots be organized into a panel? (or separated)
}
  \item{...}{
Aditional graphical arguments
}
}
\details{
Jose Luis Vicente Villardon
}
\value{
The histogram panel.
}
\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
HistogramPanel(wine[,4:7], nrows = 2, xlab="")

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
