\name{Integer2Binary}
\alias{Integer2Binary}
\title{Transforms an Integer Variable into a Binary Variable}
\description{
  Transforms an Integer Variable into a Binary Variable
}
\usage{
  Integer2Binary(y, name = "My_Factor")
}
\arguments{
  \item{y}{
Vector with the factor
}
  \item{name}{
name of the factor
}
}

\details{
Transforms an Integer vector into a Binary Indicator Matrix
}

\value{
A Binary Data Matrix
}

\author{
Jose Luis Vicente-Villardon
}

\examples{
dat=c(1, 2, 2, 4, 1, 1, 4, 2, 4)
Integer2Binary(dat,"Myfactor")
}
