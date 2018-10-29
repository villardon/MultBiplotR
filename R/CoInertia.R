Coinertia <- function(X , Y, ScalingX = 5, ScalingY = 5, dimsol=3){
  if (is.data.frame(X)) X=as.matrix(X)
  if (is.data.frame(Y)) Y=as.matrix(Y)
  
  Coiner = list()
  Coiner$Title = "Coinertia Biplot"
  Coiner$Type = "Coinertia" 
  Coiner$call <- match.call()
  Coiner$Dimension=dimsol
  #Information about X
  Coiner$XINFO$Non_Scaled_Data = X
  Coiner$XINFO$Initial_Transformation=ScalingX
  Coiner$XINFO$Means = apply(X, 2, mean)
  Coiner$XINFO$Medians = apply(X, 2, median)
  Coiner$XINFO$Deviations = apply(X, 2, sd)
  Coiner$XINFO$Minima = apply(X, 2, min)
  Coiner$XINFO$Maxima = apply(X, 2, max)
  Coiner$XINFO$P25 = apply(X, 2, quantile)[2, ]
  Coiner$XINFO$P75 = apply(X, 2, quantile)[4, ]
  Coiner$XINFO$GMean = mean(X)
  
  
  #Information about Y
  Coiner$YINFO$Non_Scaled_Data = Y
  Coiner$YINFO$Initial_Transformation=ScalingY
  Coiner$YINFO$Means = apply(Y, 2, mean)
  Coiner$YINFO$Medians = apply(Y, 2, median)
  Coiner$YINFO$Deviations = apply(Y, 2, sd)
  Coiner$YINFO$Minima = apply(Y, 2, min)
  Coiner$YINFO$Maxima = apply(Y, 2, max)
  Coiner$YINFO$P25 = apply(Y, 2, quantile)[2, ]
  Coiner$YINFO$P75 = apply(Y, 2, quantile)[4, ]
  Coiner$YINFO$GMean = mean(Y)
  
  X = InitialTransform(X, transform = ScalingX)$X
  Y = InitialTransform(X, transform = ScalingX)$X
  
  Coiner$XINFO$Scaled_Data = X
  Coiner$YINFO$Scaled_Data = Y

  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  
  # Inertia of the X matrix
  Coiner$InTotalX=sum(X^2)
  Coiner$InerRowsX=apply(X^2, 1, sum)
  Coiner$InerColsX=apply(X^2, 2, sum)
  
  # Inertia of the X matrix
  Coiner$InTotalY=sum(Y^2)
  Coiner$InerRowsY=apply(Y^2, 1, sum)
  Coiner$InerColsY=apply(Y^2, 2, sum)
  
  # Covariances among the variables of both groups
  W=(t(X) %*% Y)/I
  
  Coiner$Covariances=W
  
  # Decomposition of the covariances
  COIN=svd(W)
  U=COIN$u[,1:dimsol]
  V=COIN$v[,1:dimsol]
  D=diag(COIN$d[1:dimsol])
  Di=diag(1/COIN$d[1:dimsol])
  #Inertias
  Inertia=100*COIN$d^2 / sum(COIN$d^2)
  AcumInertia=cumsum(Inertia)
  TablaInertia=cbind(COIN$d[1:dimsol]^2, Inertia[1:dimsol], AcumInertia[1:dimsol])
  rownames(TablaInertia) = paste("Dimension", 1:dimsol)
  colnames(TablaInertia) = c("Eigenvalue", "Eplained Variance", "Cummulative")
  
  Coiner$Inertia=TablaInertia
  
  
  #Individual Biplot for X
  AX=X %*% U
  colnames(AX)=paste("Dim", 1:dimsol)
  BX=U
  rownames(BX)=colnames(X)
  colnames(BX)=paste("Dim", 1:dimsol)
  Xfitted=X %*% U %*% t(U)
  XFit=100*sum(Xfitted^2)/sum(X^2)
  
  
  #Individual Biplot for Y
  AY=Y %*% V
  BY=V
  Yfitted=Y %*% V %*% t(V)
  YFit=100*sum(Yfitted^2)/sum(Y^2)
  
  # Joint Biplot
  AXW=X %*% U %*% sqrt(Di)
  AYW=Y %*% V %*% sqrt(Di)
  BXW=U %*% sqrt(D)
  BYW=V %*% sqrt(D)
  class(Coiner)="Coinertia.SOL"
  return(Coiner)
}

CoinertiaBiplot <- function(){
  
}              

                       