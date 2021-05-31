\name{BoxPlotPanel}
\alias{BoxPlotPanel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Panel of box plots
}
\description{
Panel of box plots for a set of numerical variables and a grouping factor.
}
\usage{
BoxPlotPanel(X, groups = NULL, nrows = NULL, panel = TRUE, 
notch = FALSE, GroupsTogether = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The matrix of continuous variables
}
  \item{groups}{
The grouping factor
}
  \item{nrows}{
Number of rows of the panel.
}
  \item{panel}{
Should the plots be organized into a panel? (or separated)
}
  \item{notch}{
Should notches be used in the box plots?
}
  \item{GroupsTogether}{
Should all the groups be together in the same plot?
}
  \item{...}{
Other graphical arguments
}
}
\details{
Panel of box plots for a set of numerical variables and a grouping factor.
}
\value{
The box plot panel
}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(wine)
BoxPlotPanel(wine[,4:7], groups = wine$Origin, nrows = 2, ylab="")

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
