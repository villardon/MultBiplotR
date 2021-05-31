Biplot.PLSR1BIN <- function(plsr){
  X=plsr$X
  Y=plsr$Y
  
  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  S=dim(plsr$XScores)[2]
  
  Biplot = list()
  Biplot$Title = " PLSR - Biplot"
  Biplot$Type = "PLSR" 
  Biplot$alpha=0
  Biplot$Dimension=S
  Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$dim=S
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  if (plsr$Initial_Transformation == "Within groups standardization")  Biplot$Deviations = plsr$Deviations
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
  
  a=plsr$XScores
  b=plsr$XLoadings
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/J
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  
  Cont=CalculateContributions(plsr$ScaledX,plsr$XScores,  plsr$XLoadings )
  Biplot$Inertia=Cont$Fit*100
  Biplot$RowContributions=Cont$RowContributions
  StResponse=cor(plsr$BinaryFit$linterm, plsr$XScores)
  rownames(StResponse)="Response"
  Biplot$Structure=Cont$Structure
  Biplot$ColContributions=Cont$ColContributions
  Biplot$SupStructure=StResponse
  Biplot$SupColContributions=StResponse^2
  
  class(Biplot)="ContinuousBiplot"
  Biplot=AddBinVars2Biplot(Biplot, plsr$Y, penalization=plsr$penalization, tolerance = plsr$tolerance, maxiter = plsr$maxiter, IncludeConst=plsr$IncludeConst)
  return(Biplot)
}