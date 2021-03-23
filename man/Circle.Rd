\name{Circle}
\alias{Circle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draws a circle
}
\description{
Draws a circle for a given radius at the specified center with the given color
}
\usage{
Circle(radius = 1, origin = c(0, 0), col = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{radius}{
radius of the circle
}
  \item{origin}{
Centre of the circle
}
  \item{col}{
Color od the circle
}
  \item{...}{
Aditional graphical parameters
}
}
\details{
Draws a circle for a given radius at the specified center with the given color
}
\value{
No value is returned
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Jose Luis Vicente Villardon
}


\examples{
plot(0,0)
Circle(1,c(0,0))
}

