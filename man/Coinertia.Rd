\name{Coinertia}
\alias{Coinertia}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coinertia Analysis.
}
\description{
Calculates a Coinertia Analysis for two matrices of continuous data 
}
\usage{
Coinertia(X, Y, ScalingX = 5, ScalingY = 5, dimsol = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The first matrix in the analysis
}
  \item{Y}{
The second matrix in the analysis
}
  \item{ScalingX}{
Transformation of the X matrix
}
  \item{ScalingY}{
Transformation of the Y matrix
}
  \item{dimsol}{
Dimension of the solution
}
}
\details{
Coinertia analysis for two continuous data matrices.
}
\value{
An object of class Coinertia.SOL
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Doledec, S., & Chessel, D. (1994). Co-inertia analysis: an alternative method for studying species-environment relationships. Freshwater biology, 31(3), 277-294.

}
\author{
Jose Luis Vicente Villardon
}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

\donttest{SSI$Year == "a2006"
SSI2D=SSI[SSI$Year == "a2006",3:23]
rownames(SSI2D)=as.character(SSI$Country[SSI$Year == "a2006"])
SSIHuman2D=SSI2D[,1:9]
SSIEnvir2D=SSI2D[,10:16]
SSIEcon2D=SSI2D[,17:21]
Coin=Coinertia(SSIHuman2D, SSIEnvir2D)}

}
