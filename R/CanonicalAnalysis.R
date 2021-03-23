Canonical.Variate.Analysis <- function(X, group, InitialTransform = 5){
  res=CanonicalBiplot(X, group=group, InitialTransform = InitialTransform, LDA=TRUE, MANOVA = TRUE)
  unclass(res)
  res$Title="Canonical Analysis"
  class(res)="CVA"
  return(res)
}

plot.CVA<- function(x, A1=1, A2=2, ...){
  plot.Canonical.Biplot(x, A1=A1, A2=A2, PlotVars=FALSE, ...)
}

summary.CVA<-function(object, ...){
  unclass(object)
  class(object)="Canonical.Biplot"
  summary(object)
}
 