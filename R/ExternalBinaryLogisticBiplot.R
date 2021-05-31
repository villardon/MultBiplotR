
logit <- function(p) {
  logit = log(p/(1 - p))
  return(logit)
}

ExternalBinaryLogisticBiplot <- function(Pco, IncludeConst=TRUE, penalization=0.2, freq=NULL, tolerance = 1e-05, maxiter = 100) {
  if (!(Pco$TypeData=="Binary")) stop("Data must be Binary for a External Binary Logistic Biplot")
  n=dim(Pco$Data)[1]
  p=dim(Pco$Data)[2]
  dimens=dim(Pco$RowCoordinates)[2]
  x=Pco$RowCoordinates

  Pco$ColumnParameters=matrix(0,p,dimens+1)
  Res=list()
  Res$Deviances=matrix(0,p,1)
  Res$Dfs=matrix(0,p,1)
  Res$pvalues=matrix(0,p,1)
  Res$Bonferroni=matrix(0,p,1)
  Res$Nagelkerke=matrix(0,p,1)
  Res$R2=matrix(0,p,1)
  Res$PercentsCorrec=matrix(0,p,1)
  Pco$DevianceTotal=0
  Pco$p=1
  Pco$TotalPercent=0
  
  for (i in 1:p){
    y=Pco$Data[,i]
    fit=RidgeBinaryLogistic(y,x,tolerance = tolerance, maxiter = maxiter, penalization=penalization, cte=IncludeConst)
    Pco$ColumnParameters[i,]=t(fit$beta)
    Res$Deviances[i]=fit$Dif
    Res$Dfs[i]=fit$df
    Res$pvalues[i]=fit$p
    Res$Bonferroni[i]=(fit$p * p)* ((fit$p * p)<=1) + (((fit$p * p)>1))
    Res$Nagelkerke[i]=fit$Nagelkerke
    Res$R2[i]=fit$R2
    Res$PercentsCorrec[i]=fit$PercentCorrect
    Pco$TotalPercent=Pco$TotalPercent+sum(y==fit$Prediction)
  }
  rownames(Pco$ColumnParameters)=colnames(Pco$Data)
  colnames(Pco$ColumnParameters)=paste("b",0:dimens, sep="")
  
  esp = cbind(rep(1,n), Pco$RowCoordinates) %*% t(Pco$ColumnParameters)
  pred = exp(esp)/(1 + exp(esp))
  
  esp0 = matrix(rep(1,n), n,1) %*% Pco$ColumnParameters[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))
  
  d1 = -2 * apply(Pco$Data * log(pred0) + (1 - Pco$Data) * log(1 - pred0),2,sum)
  d2 = -2 * apply(Pco$Data * log(pred) + (1 - Pco$Data) * log(1 - pred),2,sum)
  
  d = d1 - d2
  
  dd = sqrt(rowSums(cbind(1,Pco$ColumnParameters[, 2:(dimens + 1)])^2))
  Res$Loadings = diag(1/dd) %*% Pco$ColumnParameters[, 2:(dimens + 1)]
  Res$Tresholds = Pco$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)
  
  Pco$TotalPercent=Pco$TotalPercent/(n*p)
  Pco$DevianceTotal=sum(Res$Deviances)
  Pco$TotalDf=sum(Res$Dfs)
  Pco$p=1-pchisq(Pco$DevianceTotal, df = Pco$TotalDf)
  Res=as.data.frame(Res)
  rownames(Res)=colnames(Pco$Data)
  Pco$VarInfo=Res
  Pco$ClusterType="us"
  Pco$Clusters = as.factor(matrix(1,nrow(Pco$RowCoordinates), 1))
  Pco$ClusterColors="blue"
  Pco$ClusterNames="Cluster1"
  class(Pco)="External.Binary.Logistic.Biplot"
  return(Pco)
}
