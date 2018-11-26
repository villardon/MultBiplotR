plot.PCA.Analysis <- function(x, A1 = 1, A2 = 2, CorrelationCircle=FALSE, ...){
  plot.ContinuousBiplot(x, A1=A1, A2=A2, PlotVars=FALSE, ...)
  if (CorrelationCircle) CorrelationCircle(x) 
  else plot.ContinuousBiplot(x, A1=A1, A2=A2, PlotVars=FALSE, ...)
}