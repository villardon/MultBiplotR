\name{InBox}
\alias{InBox}

\title{
Checks if a point is inside a box.
}
\description{
Checks if a point is inside a box. The point is specified bi its x and y coordinates and the bom with the minimum and maximum values on both coordinate axis: xmin, xmax, ymin, ymax. The vertices of the box are then (xmin, ymin), (xmax, ymin), (xmax, ymax) and (xmin, ymax)
}
\usage{
InBox(x, y, xmin, xmax, ymin, ymax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
x coordinate of the point
}
  \item{y}{
x coordinate of the point
}
  \item{xmin}{
minimum value of X
}
  \item{xmax}{
maximum value of X
}
  \item{ymin}{
minimum value of Y
}
  \item{ymax}{
maximum value of Y
}
}

\value{
Returns a logical value : TRUE if the point is inside the box and FALSE otherwise.
}

\author{
Jose Luis Vicente Villardon
}



\examples{
InBox(0, 0, -1, 1, -1, 1) 
}