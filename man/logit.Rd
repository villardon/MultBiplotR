\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logit function
}
\description{
Calculates the logit of a probability
}
\usage{
logit(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
A probability
}
}
\details{
Calculates the logit of a probability
}
\value{
The lo git of the provided probebility
}

\author{
Jose Luis Vicente Villardón
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p) 
{
    logit = log(p/(1 - p))
    return(logit)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
