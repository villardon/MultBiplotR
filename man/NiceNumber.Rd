\name{NiceNumber}
\alias{NiceNumber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nice numbers: simple decimal numbers
}
\description{
Calculates a close nice number, i. e. a number with simple decimals.
}
\usage{
NiceNumber(x = 6, round = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A number
}
  \item{round}{
Should the number be rounded?
}
}
\details{
Calculates a close nice number, i. e. a number with simple decimals.
}
\value{
A number with simple decimals
}
\references{
Heckbert, P. S. (1990). Nice numbers for graph labels. In Graphics Gems (pp. 61-63). Academic Press Professional, Inc..
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
\code{\link{PrettyTicks}}
}
\examples{
NiceNumber(0.892345)
}