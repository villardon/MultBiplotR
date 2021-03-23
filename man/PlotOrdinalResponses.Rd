\name{PlotOrdinalResponses}
\alias{PlotOrdinalResponses}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the response functions along the directions of best fit.
}
\description{
Plot the response functions along the directions of best fit for the selected dimensions
}
\usage{
PlotOrdinalResponses(olb, A1 = 1, A2 = 2, inf = -12, sup = 12, 
Legend = TRUE, WhatVars=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{olb}{
An object of class "Ordinal.Logistic.Biplot"}
  \item{A1}{
First dimension of the plot.
}
  \item{A2}{
Second dimension of the plot
}
  \item{inf}{
Lower limit of the representation
}
  \item{sup}{
Upper limit of the representation
}
  \item{Legend}{
Should a legend be plotted
}
  \item{WhatVars}{
A vector with the numbers of the variables to be plotted. If NULL all the variables are plotted.
}
}
\details{
Plot the response functions along the directions of best fit for the selected dimensions
}
\value{
A plot describing the behaviour of the variable 
}

\author{
Jose Luis Vicente Villardon
}
\examples{

data(Doctors)
    olb = OrdLogBipEM(Doctors,dim = 2, nnodes = 10, initial=4,  tol = 0.001, 
    maxiter = 100, penalization = 0.1, show=TRUE)
    PlotOrdinalResponses(olb, WhatVars=c(1,2,3,4))
}
