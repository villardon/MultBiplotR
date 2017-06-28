\name{plot.CCA.sol}
\alias{plot.CCA.sol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the solution of a Canonical Correspondence Analysisis
}
\description{
Plots the solution of a Canonical Correspondence Analysisis using similar parameters to the continuous biplot
}
\usage{
\method{plot}{CCA.sol}(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0,
                 PlotSites = TRUE, PlotSpecies = TRUE, PlotEnv = TRUE,
                 LabelSites = TRUE, LabelSpecies = TRUE, LabelEnv =
                 TRUE, TypeSites = "wa", SpeciesQuality = FALSE,
                 MinQualityVars = 0.3, dp = 0, pr = 0, PlotAxis =
                 FALSE, TypeScale = "Complete", ValuesScale =
                 "Original", mode = "a", CexSites = NULL, CexSpecies =
                 NULL, CexVar = NULL, ColorSites = NULL, ColorSpecies =
                 NULL, ColorVar = NULL, PchSites = NULL, PchSpecies =
                 NULL, PchVar = NULL, SizeQualSites = FALSE,
                 SizeQualSpecies = FALSE, SizeQualVars = FALSE,
                 ColorQualSites = FALSE, ColorQualSpecies = FALSE,
                 ColorQualVars = FALSE, SmartLabels = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{CCA.sol} here~~
}
  \item{A1}{
%%     ~~Describe \code{A1} here~~
}
  \item{A2}{
%%     ~~Describe \code{A2} here~~
}
  \item{ShowAxis}{
%%     ~~Describe \code{ShowAxis} here~~
}
  \item{margin}{
%%     ~~Describe \code{margin} here~~
}
  \item{PlotSites}{
%%     ~~Describe \code{PlotSites} here~~
}
  \item{PlotSpecies}{
%%     ~~Describe \code{PlotSpecies} here~~
}
  \item{PlotEnv}{
%%     ~~Describe \code{PlotEnv} here~~
}
  \item{LabelSites}{
%%     ~~Describe \code{LabelSites} here~~
}
  \item{LabelSpecies}{
%%     ~~Describe \code{LabelSpecies} here~~
}
  \item{LabelEnv}{
%%     ~~Describe \code{LabelEnv} here~~
}
  \item{TypeSites}{
%%     ~~Describe \code{TypeSites} here~~
}
  \item{SpeciesQuality}{
%%     ~~Describe \code{SpeciesQuality} here~~
}
  \item{MinQualityVars}{
%%     ~~Describe \code{MinQualityVars} here~~
}
  \item{dp}{
%%     ~~Describe \code{dp} here~~
}
  \item{pr}{
%%     ~~Describe \code{pr} here~~
}
  \item{PlotAxis}{
%%     ~~Describe \code{PlotAxis} here~~
}
  \item{TypeScale}{
%%     ~~Describe \code{TypeScale} here~~
}
  \item{ValuesScale}{
%%     ~~Describe \code{ValuesScale} here~~
}
  \item{mode}{
%%     ~~Describe \code{mode} here~~
}
  \item{CexSites}{
%%     ~~Describe \code{CexSites} here~~
}
  \item{CexSpecies}{
%%     ~~Describe \code{CexSpecies} here~~
}
  \item{CexVar}{
%%     ~~Describe \code{CexVar} here~~
}
  \item{ColorSites}{
%%     ~~Describe \code{ColorSites} here~~
}
  \item{ColorSpecies}{
%%     ~~Describe \code{ColorSpecies} here~~
}
  \item{ColorVar}{
%%     ~~Describe \code{ColorVar} here~~
}
  \item{PchSites}{
%%     ~~Describe \code{PchSites} here~~
}
  \item{PchSpecies}{
%%     ~~Describe \code{PchSpecies} here~~
}
  \item{PchVar}{
%%     ~~Describe \code{PchVar} here~~
}
  \item{SizeQualSites}{
%%     ~~Describe \code{SizeQualSites} here~~
}
  \item{SizeQualSpecies}{
%%     ~~Describe \code{SizeQualSpecies} here~~
}
  \item{SizeQualVars}{
%%     ~~Describe \code{SizeQualVars} here~~
}
  \item{ColorQualSites}{
%%     ~~Describe \code{ColorQualSites} here~~
}
  \item{ColorQualSpecies}{
%%     ~~Describe \code{ColorQualSpecies} here~~
}
  \item{ColorQualVars}{
%%     ~~Describe \code{ColorQualVars} here~~
}
  \item{SmartLabels}{
%%     ~~Describe \code{SmartLabels} here~~
}
  \item{...}{
Aditional graphical parameters.
}
}
\details{
The plotting procedure is similar to the one used for continuous biplots including the calibration of the environmental variables.
}
\value{
No value returned
}
\references{
CCA
}
\author{
Jose Luis Vicente Villardon
}

\seealso{
 \code{\link{plot.ContinuousBiplot}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
}

