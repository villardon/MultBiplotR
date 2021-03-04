ConstrainedLogisticBiplot <- function(Y, X, dim=2, Scaling=5,  tolerance = 1e-05, maxiter = 100, penalization = 0.1){
  Y=as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  r=min(p,q)
  X=as.matrix(X)
  Xo=X
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(Scaling)) 
    Scaling = ContinuousDataTransform[Scaling]
  

  Data = InitialTransform(X, transform = Scaling)
  X = Data$X
  if (Scaling=="Within groups standardization") Biplot$Deviations = Data$ColStdDevs
  
  Yh=matrix(0, n,p)
  Betas=matrix(0, q+1, p)
  X=cbind(rep(1,n),X)
  for (i in 1:p){
    ajuste=RidgeBinaryLogistic(Y[,i], X, tolerance = tolerance, maxiter = maxiter, penalization = penalization)
    Yh[,i]=ajuste$linterm
    Betas[,i]=ajuste$beta
  }
  Betas=t(Betas)
  rownames(Betas)=colnames(Y)
  colnames(Betas)=colnames(X)
  colnames(Betas)[1]="Intercept"
  Yh = InitialTransform(Yh, transform = 5)$X
  SD = svd(Yh)
  EV = SD$d[1:r]^2
  Inertia = round((EV/sum(EV)) * 100, digits = 3)
  CumInertia = cumsum(Inertia)
  a = SD$u[,1:dim] %*% diag(SD$d[1:dim])
  rownames(a)=rownames(Y)
  freq = matrix(1, n, 1)
  aa=cbind(matrix(1, n, 1), a)
  b = matrix(0, p, dim + 1)
  for (j in 1:p) {
    b[j, ] = RidgeBinaryLogisticFit(Y[, j], aa, freq = matrix(1, n, 1), tolerance=tolerance, maxiter=maxiter)
  }
  
  rownames(b)=colnames(Y)
  colnames(b)=c("intercept", paste("dim",1:dim))
  
  esp = aa %*% t(b)
  pred = exp(esp)/(1 + exp(esp))
  loglikelyhood = sum(sum(Y * log(pred) + (1 - Y) * log(1 - pred)))
  pred2 = matrix(as.numeric(pred > 0.5), n, p)
  acier = matrix(as.numeric(Y == pred2), n, p)
  acierfil = 100*apply(acier,1,sum)/p
  aciercol = 100*apply(acier,2,sum)/n
  
  presences=apply(Y, 2, sum)
  absences=n-presences
  sens = apply((acier==1) & (Y==1), 2, sum)/presences
  spec = apply((acier==1) & (Y==0), 2, sum)/absences
  totsens = sum((acier==1) & (Y==1))/sum(presences)
  totspec = sum((acier==1) & (Y==0))/sum(absences)

  
  gfit = (sum(sum(acier))/(n * p)) * 100
  
  esp0 = freq %*% b[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))
  
  d1 = -2 * apply(Y * log(pred0) + (1 - Y) * log(1 - pred0),2,sum)
  d2 = -2 * apply(Y * log(pred) + (1 - Y) * log(1 - pred),2,sum)
  
  d = d1 - d2
  ps = matrix(0, p, 1)
  for (j in 1:p) ps[j] = 1 - pchisq(d[j], 1)
  
  resid = Y - pred
  r2 = 1 - (sum(sum(resid^2))/sum(sum(Y^2)))
  pred2 = (esp > 0.5)
  corrfil = acierfil
  corrcol = aciercol
  rownames(b)=colnames(Y)
  
  Res=list()
  Res$Biplot="Binary Logistic (Constrained)"
  Res$Type= "Binary Logistic (Constrained)"
  Res$Betas=Betas
  Res$RowCoordinates=a
  Res$ColumnParameters=b
  
  Res$NullDeviances=d1
  Res$ModelDeviances=d2
  Res$Deviances=d
  Res$Dfs=rep(dim, p)
  Res$pvalues=ps
  Res$CoxSnell=1-exp(-1*Res$Deviances/n)
  Res$Nagelkerke=Res$CoxSnell/(1-exp((Res$NullDeviances/(-2)))^(2/n))
  Res$MacFaden=1-(Res$ModelDeviances/Res$NullDeviances)
  
  Res$TotalPercent=gfit
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$DevianceTotal=sum(Res$Deviances)
  
  dd = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(dim + 1)])^2))
  Res$Loadings = solve(diag(dd)) %*% Res$ColumnParameters[, 2:(dim + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)
  
  nn=n*p
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  
  Res$R2 = apply((Y-esp)^2,2, sum)/apply((Y)^2,2, sum)

  Res$TotR2 = sum((Y-esp)^2) /sum((Y)^2)
  pred= matrix(as.numeric(esp>0.5),n , p)
  verdad = matrix(as.numeric(Y==pred),n , p)
  Res$PercentsCorrec=apply(verdad, 2, sum)/n
  Res$TotalPercent=sum(verdad)/(n*p)
  Res$Sensitivity=sens
  Res$Specificity=spec
  Res$TotalSensitivity=totsens
  Res$TotalSpecificity=totspec
  Res$TotalDf = dim*p
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)
  
  Res$ClusterType="us"
  Res$Clusters = as.factor(matrix(1,n, 1))
  Res$ClusterColors="blue"
  Res$ClusterNames="ClusterTotal"
  class(Res) = "Binary.Logistic.Biplot"
  Res=AddContVars2Biplot(Res,  Xo, dims=dim, Scaling = Scaling)
  return(Res)
}

