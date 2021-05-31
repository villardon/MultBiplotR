\name{SSI}
\alias{SSI}
\docType{data}
\title{
Sustainability Society Index
}
\description{
Sustainability Society Index
}
\usage{data("SSI")}
\format{
  A data frame with 924 observations on the following 23 variables.
  \describe{
    \item{\code{Year}}{a factor with levels \code{a2006} \code{a2008} \code{a2010} \code{a2012} \code{a2014} \code{a2016}}

    
    \item{\code{Country}}{a factor with levels \code{Albania} \code{Algeria} 
    \code{Angola} \code{Argentina} \code{Armenia} \code{Australia} \code{Austria}
    \code{Azerbaijan} \code{Bangladesh} \code{Belarus} \code{Belgium} \code{Benin} 
    \code{Bhutan} \code{Bolivia} \code{Bosnia-Herzegovina} \code{Botswana} \code{Brazil}
    \code{Bulgaria} \code{Burkina_Faso} \code{Burundi} \code{Cambodia} \code{Cameroon} 
    \code{Canada} \code{Central_African_Republic} \code{Chad} \code{Chile} 
    \code{China} \code{Colombia} \code{Congo} \code{Congo_Democratic_Rep.} 
    \code{Costa_Rica} \code{Cote_dIvoire} \code{Croatia} \code{Cuba}
    \code{Cyprus} \code{Czech_Republic} \code{Denmark} \code{Dominican_Republic} 
    \code{Ecuador} \code{Egypt} \code{El_Salvador} \code{Estonia} \code{Ethiopia} 
    \code{Finland} \code{France} \code{Gabon} \code{Gambia} \code{Georgia} 
    \code{Germany} \code{Ghana} \code{Greece} \code{Guatemala} \code{Guinea} 
    \code{Guinea-Bissau} \code{Guyana} \code{Haiti} \code{Honduras} 
    \code{Hungary} \code{Iceland} \code{India} \code{Indonesia} \code{Iran} 
    \code{Iraq} \code{Ireland} \code{Israel} \code{Italy} \code{Jamaica} 
    \code{Japan} \code{Jordan} \code{Kazakhstan} \code{Kenya} \code{Korea._North} 
    \code{Korea._South} \code{Kuwait} \code{Kyrgyz_Republic} \code{Laos} 
    \code{Latvia} \code{Lebanon} \code{Lesotho} \code{Liberia} \code{Libya} 
    \code{Lithuania} \code{Luxembourg} \code{Macedonia} \code{Madagascar} 
    \code{Malawi} \code{Malaysia} \code{Mali} \code{Malta} \code{Mauritania} 
    \code{Mauritius} \code{Mexico} \code{Moldova} \code{Mongolia} \code{Montenegro} 
    \code{Morocco} \code{Mozambique} \code{Myanmar} \code{Namibia} \code{Nepal} 
    \code{Netherlands} \code{New_Zealand} \code{Nicaragua} \code{Niger} \code{Nigeria} 
    \code{Norway} \code{Oman} \code{Pakistan} \code{Panama} \code{Papua_New_Guinea} 
    \code{Paraguay} \code{Peru} \code{Philippines} \code{Poland} \code{Portugal} 
    \code{Qatar} \code{Romania} \code{Russia} \code{Rwanda} \code{Saudi_Arabia} 
    \code{Senegal} \code{Serbia} \code{Sierra_Leone} \code{Singapore} 
    \code{Slovak_Republic} \code{Slovenia} \code{South_Africa} \code{Spain} 
    \code{Sri_Lanka} \code{Sudan} \code{Sweden} \code{Switzerland} \code{Syria} 
    \code{Taiwan} \code{Tajikistan} \code{Tanzania} \code{Thailand} \code{Togo} \code{Trinidad_and_Tobago} \code{Tunisia} \code{Turkey} \code{Turkmenistan} 
    \code{Uganda} \code{Ukraine} \code{United_Arab_Emirates} \code{United_Kingdom} \code{United_States} \code{Uruguay} \code{Uzbekistan} \code{Venezuela} \code{Vietnam} \code{Yemen} \code{Zambia} \code{Zimbabwe}}
    
    \item{\code{Sufficient_Food}}{a numeric vector}
    \item{\code{Sufficient_to_Drink}}{a numeric vector}
    \item{\code{Safe_Sanitation}}{a numeric vector}
    \item{\code{Education_}}{a numeric vector}
    \item{\code{Healthy_Life}}{a numeric vector}
    \item{\code{Gender_Equality}}{a numeric vector}
    \item{\code{Income_Distribution}}{a numeric vector}
    \item{\code{Population_Growth}}{a numeric vector}
    \item{\code{Good_Governance}}{a numeric vector}
    \item{\code{Biodiversity_}}{a numeric vector}
    \item{\code{Renewable_Water_Resources}}{a numeric vector}
    \item{\code{Consumption}}{a numeric vector}
    \item{\code{Energy_Use}}{a numeric vector}
    \item{\code{Energy_Savings}}{a numeric vector}
    \item{\code{Greenhouse_Gases}}{a numeric vector}
    \item{\code{Renewable_Energy}}{a numeric vector}
    \item{\code{Organic_Farming}}{a numeric vector}
    \item{\code{Genuine_Savings}}{a numeric vector}
    \item{\code{GDP}}{a numeric vector}
    \item{\code{Employment}}{a numeric vector}
    \item{\code{Public_Debt}}{a numeric vector}
  }
}
\details{
Sustainability Society Index
}
\source{
https://ssi.wi.th-koeln.de
}
\references{
Gallego-Alvarez, I., Galindo-Villardon, M. P., & Rodriguez-Rosa, M. (2015). Analysis of the Sustainable Society Index Worldwide: A Study from the Biplot Perspective. Social Indicators Research, 120(1), 29-65. https://doi.org/10.1007/s11205-014-0579-9
}
\examples{
data(SSI)
## maybe str(SSI) ; plot(SSI) ...
}
\keyword{datasets}
