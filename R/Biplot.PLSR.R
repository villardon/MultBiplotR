Biplot.PLSR <- function(plsr){
  X=plsr$X
  Y=plsr$Y
  
  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  S=dim(plsr$XScores)[2]
  
  Biplot = list()
  Biplot$Title = " PLSR - Biplot"
  Biplot$Type = "PLSR" 
  Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$Dimension=S
  Biplot$alpha=0
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
  
  a = plsr$XScores
  b = plsr$XLoadings
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/J
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  Cont=CalculateContributions(plsr$ScaledX,plsr$XScores,  plsr$XLoadings)
  Biplot$EigenValues=Cont$Fit
  Biplot$Inertia=Cont$Fit*100
  Biplot$CumInertia=cumsum(Biplot$Inertia)
  Biplot$RowContributions=Cont$RowContributions
  Biplot$ColContributions=Cont$ColContributions
  Biplot$Structure=Cont$Structure
  class(Biplot)="ContinuousBiplot"
  
  # Adding Y information as a supplementary Biplot
  
  YBiplot=list()
  YBiplot$Title = " PLSR - Biplot (Y)"
  YBiplot$Type = "PLSR" 
  YBiplot$Initial_Transformation=plsr$Initial_Transformation
  YBiplot$ncols=K
  YBiplot$Means = apply(Y, 2, mean)
  YBiplot$Medians = apply(Y, 2, median)
  YBiplot$Deviations = apply(Y, 2, sd)
  YBiplot$Minima = apply(Y, 2, min)
  YBiplot$Maxima = apply(Y, 2, max)
  YBiplot$P25 = apply(Y, 2, quantile)[2, ]
  YBiplot$P75 = apply(Y, 2, quantile)[4, ]
  YBiplot$b0 = rep(0,K)
  
  YBiplot$ColCoordinates = plsr$YWeights/scf
  if (K>1)
  Cont=CalculateContributions(plsr$ScaledY,plsr$YScores,  plsr$YLoadings)
  YBiplot$ColContributions=Cont$ColContributions
  YBiplot$Structure=Cont$Structure
  class(YBiplot)="ContSupVarsBiplot"
  Biplot$ContSupVarsBiplot = YBiplot
  return(Biplot)
}
