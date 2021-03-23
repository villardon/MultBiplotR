\name{Matrix2Proximities}
\alias{Matrix2Proximities}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix to Proximities
}
\description{
Converts a matrix of proximities into a Proximities object as used in Principal Coordinates or MDS
}
\usage{
Matrix2Proximities(x, TypeData = "User Provided", 
Type = c("dissimilarity", "similarity", "products"), 
Coefficient = "None", Transformation = "None", Data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The matrix of proximities (a symmetrical matrix)
}
  \item{TypeData}{
By default is User provided but could be any type.
}
  \item{Type}{
Type of proximity: dissimilarity, similarity or scalar product. If not provided, the default is dissimilarity
}
  \item{Coefficient}{
Name of the procedure to calculate the proximities (if any).
}
  \item{Transformation}{
Transformation used to calculate dissimilarities from similarities (if any)
}
  \item{Data}{
Raw data used to calculate the proximity (if any).
}
}
\details{
Converts a matrix of proximities into a Proximities object as used in Principal Coordinates or MDS aading some extra information about the procedure used to obtain the proximities. Is mainly used when the proximities matrix has been provided by the user and not calculated from raw data using BinaryProximities, ContinuousDistances or any other function.
}
\value{
An object of class \code{Proximities} containing the proximities matrix and some extra information about it.
}

\author{
Jose Luis Vicente Villardon
}

