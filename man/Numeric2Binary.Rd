\name{Numeric2Binary}
\alias{Numeric2Binary}

\title{Converts a numeric variable into a binary one}
\description{
Converts a numeric variable into a binary one using a cut point
}
\usage{
Numeric2Binary(y, name= "MyVar", cut = NULL)
}

\arguments{
  \item{y}{
Vector containing the numeric values
}
  \item{name}{
Name of the variable
}
  \item{cut}{
Cut point to cut the values of the variable. If is NULL the median is used.
}
}
\details{
Converts a numeric variable into a binary one using a cut point. If the cut is NULL the median is used.
}
\value{
A binary Variable
}

\author{
Jose Luis Vicente-Villardon
}

\seealso{
Dataframe2BinaryMatrix
}
\examples{
y=c(1, 1.2, 3.2, 2.4, 1.7, 2.2, 2.7, 3.1)
Numeric2Binary(y)
}

