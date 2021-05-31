\name{Games_Howell}
\alias{Games_Howell}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Games-Howell post-hoc tests for Welch's one-way analysis
}
\description{
This function produces results from Games-Howell post-hoc tests for Welch's one-way analysis of variance (ANOVA) for a matrix of numeric data and a grouping variable.
}
\usage{
Games_Howell(data, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The matrix of continuous data.
}
  \item{group}{
The grouping variable
}
}
\details{
This function produces results from Games-Howell post-hoc tests for Welch's one-way analysis of variance (ANOVA) for a matrix of numeric data and a grouping variable.
}
\value{
The tests for each column of the data matrix
}
\references{
Ruxton, G. D., & Beauchamp, G. (2008). Time for some a priori thinking about post hoc testing. Behavioral ecology, 19(3), 690-693.
}
\author{
Jose Luis Vicente Villardon
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
