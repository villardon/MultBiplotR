\encoding{UTF-8}
\name{NominalDistances}
\alias{NominalDistances}

\title{
Distances among individuals with nominal variables 
}
\description{
This function computes several measures of distance (or similarity) among individuals from a nominal data matrix. 
}
\usage{
NominalDistances(X, method = 1, diag = FALSE, upper = FALSE, similarity = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Matrix or data.frame with the nominal variables.
}
  \item{method}{
An integer between 1 and 6. See details 
}
  \item{diag}{
A logical value indicating whether the diagonal of the distance matrix should be printed.
}
  \item{upper}{
a logical value indicating whether the upper triangle of the distance matrix should be printed.
}
  \item{similarity}{
A logical value indicating whether the similarity matrix should be computed.
}
}
\details{
Let be the table of nominal data. All these distances are of type \eqn{d=\sqrt{1-s}}{d = sqrt(1 - s)} with \emph{s} a similarity coefficient.
\describe{

\item{1 = Overlap method}{The overlap measure simply counts the number of attributes that match in the two data instances.}
\item{2 = Eskin}{Eskin et al. proposed a normalization kernel for record-based network intrusion detection data. The original measure is distance-based and assigns a weight of \eqn{\frac{2}{n_{k}^{2}}} for mismatches; when adapted to similarity, this becomes a weight of \eqn{\frac{n_{k}^{2}}{n_{k}^{2}+2}}.This measure gives more weight to mismatches that occur on attributes that take many values.}
\item{3=IOF (Inverse Occurrence Frequency .)}{
This measure assigns lower similarity to mismatches on more frequent values. The IOF measure is related to the concept of inverse document frequency which comes from information retrieval, where it is used to signify the relative number of documents that contain a spe- cific word.}
\item{4 = OF (Ocurrence Frequency)}{This measure gives the opposite weighting of the IOF measure for mismatches, i.e., mismatches on less frequent values are assigned lower similarity and mismatches on more frequent values are assigned higher similarity}
\item{5 = Goodall3}{This measure assigns a high similarity if the matching values are infrequent regardless of the frequencies of the other values.}
\item{6 = Lin}{This measure gives higher weight to matches on frequent values, and lower weight to mismatches on infrequent values.}
           
}
}
\value{
An object of class distance
}
\references{
Boriah, S., Chandola, V. & Kumar,V.(2008). Similarity measures for categorical data: A comparative evaluation. In proceedings of the eight SIAM International Conference on Data Mining, pp 243--254.

}
\author{
Jose L. Vicente-Villardon
}


\seealso{
\code{\link{BinaryDistances}},\code{\link{ContinuousDistances}}
}
\examples{
\dontrun{
data(Env)
Distance<-NominalDistances(Env,upper=TRUE,diag=TRUE,similarity=FALSE,method=1)
}
}

\keyword{Nominal}
\keyword{Distance}% __ONLY ONE__ keyword per line
\keyword{Multivariate}