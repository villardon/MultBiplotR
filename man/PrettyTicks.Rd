\name{PrettyTicks}
\alias{PrettyTicks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates loose axis ticks and labels using nice numbers
}
\description{
Calculates axis ticks and labels using nice numbers
}
\usage{
PrettyTicks(min = -3, max = 3, ntick = 5)
}

\arguments{
  \item{min}{
Minimum value on the axis
}
  \item{max}{
maximum value on the axis.
}
  \item{ntick}{
Approximated number of desired ticks
}
}
\details{
Calculates axis ticks and labels using nice numbers. The resulting labels are known as loose labels.
}
\value{
A list with the following fields
 \item{ticks }{Ticks for the axis}
\item{labels }{The corresponding labels}
}
\references{
Heckbert, P. S. (1990). Nice numbers for graph labels. In Graphics Gems (pp. 61-63). Academic Press Professional, Inc..
}
\author{
Jose Luis Vicente Villardon
}
\seealso{
\code{\link{NiceNumber}}
}
\examples{
PrettyTicks(-4, 4, 5)
}

