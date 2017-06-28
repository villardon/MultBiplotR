\name{Doctors}
\alias{Doctors}
\docType{data}
\title{
Data set extracted from the Careers of doctorate holders survey carried out by Spanish Statistical Office in 2008.
}
\description{
The sample data, as part of a large survey, corresponds to 100 people who have the PhD degree
and it shows the level of satisfaction of the doctorate holders about some issues.
}
\usage{data(Doctors)}
\format{
  This data frame contains 100 observation for the following 5 ordinal variables,
  with four categories each:   (1= "Very Satisfied", 2= "Somewhat Satisfied",3="Somewhat dissatisfied",  4="Very dissatisfied")
  \describe{
    \item{Salary}{}
    \item{Benefits}{}
    \item{Job Security}{}
    \item{Job Location}{}
    \item{Working conditions}{}
  }
}
\source{
Spanish Statistical Institute. Survey of PDH holders, 2006. 
URL: http://www.ine.es.
}
\examples{
data(Doctors)
## maybe str(Doctors) ; plot(Doctors) ...
}
\keyword{datasets}
