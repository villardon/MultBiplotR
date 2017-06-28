\name{SeparateVarTypes}
\alias{SeparateVarTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Separation of different types of variables into a list
}
\description{
The procedure creates a list in which each field contains the variables of the same type.
}
\usage{
SeparateVarTypes(X, TypeVar = NULL, TypeFit = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A data frame
}
  \item{TypeVar}{
A vector of characters defining the type of each variable. If not provided the procedure tries to gess the type of each variable. See details for types
}
  \item{TypeFit}{
A vector of characters defining the type of fit for each variable. If not provided the procedure tries to gess the type of fit for each variable. See details for types
}
}
\details{
The procedure creates a list in which each field contains the variables of the same type.
The type of Variable can be specified in a vector TypeVar and the type of fit in a vector TypeFit.
The TypeVar is a vector of characters with as many components as variables with types coded as:

"c" - Continuous (1)

"b" - Binary (2)

"n" - Nominal (3)

"o" - Ordinal (4)

"f" - Frequency (5)

"a" - Abundance (5)

Numbers rhather than characters can also be used.
Unless specified in TypeVar, numerical variables are "Continuous", factors are "Nominal", ordered factors are "Ordinal". Factors with just two values 
are considered as "Binary".
"Frequencies" and "abundances" should be specified by the user. If Typevar has length 1, all the variables are supposed to have the same type.

The typeFit is a vector of characters containing the type of fit used for each variab, coded as:

"a" - Average (1)

"wa" - Weighted Average (2)

"r" - Regression (Linear or logistic depending on the type of variable) (3)

"g" - Gaussian (Equal tolerances) (4)

"g1" - Gaussian (Different tolerances)  (5)


Numbers rhather than characters can also be used. Unless specified numerical variables are fitted with linear regression, 
factors with logistic biplots, frequencies with weighted averages and abundances with gaussian regression. 

}
\value{
A list containing the following fields
\item{Continuous }{A list containing a data frame with the numeric variables and a character vector with the type of fit for each variable}
\item{Binary }{A list containing a data frame with the binary variables and a character vector with the type of fit for each variable}
\item{Nominal }{A list containing a data frame with the nominal variables and a character vector with the type of fit for each variable}
\item{Ordinal }{A list containing a data frame with the ordinal variables and a character vector with the type of fit for each variable}
\item{Frequency }{A list containing a data frame with the frequency variables and a character vector with the type of fit for each variable}
\item{Abundance }{A list containing a data frame with the abundance variables and a character vector with the type of fit for each variable}

}

\author{
Jose Luis Vicente Villardon
}

\examples{
data(Protein)
SepData=SeparateVarTypes(Protein)
SepData
}

