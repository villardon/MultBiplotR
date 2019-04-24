BinaryLogBiplotEM <- function(x, freq=matrix(1,nrow(x),1), aini=NULL, dimens = 2, nnodos = 15, tol = 1e-04, maxiter = 100, penalization = 0.2) {
  # Nodos y ponderaciones de la cuadratura de Gauss-Hermite
  Q = Multiquad(nnodos, dimens)
  X = Q$X
  A = Q$A
  q = dim(X)[1]
  s = dim(x)[1]
  
  if (is.null(freq)) freq=matrix(1,s,1)
  n = dim(x)[2]
  # Parametros iniciales
  print("Calculating initial parameters")
  par = matrix(0, n, dim + 1)
  
  xx=as.matrix(cbind(x, 1 - x))
  # corr=CA(xx,dim=dim)
  # ability=corr$RowCoordinates[,1:dim]
  
  if (is.null(aini)){
    mod2 <- mirt(x,model=dim)
    ability <- fscores(mod2,method = "EAP", full.scores = TRUE)}
  else
    ability=aini
  
  print("  ")
  rownames(ability)=rownames(x)
  colnames(ability)=paste("Axis",1:dim)
  
  for (j in 1:n) par[j, ] = RidgeBinaryLogisticFit(x[, j], cbind(matrix(1, s, 1), ability), freq, tol, maxiter, penalization)
  lin = cbind(matrix(1, s, 1), ability) %*% t(par)
  PI = 1/(1 + exp(-1 * lin))
  logLikold = sum(diagonal(freq) %*% log(PI^x * ((1 - PI)^(1 - x))))
  
  # Inicializacion de los parámetros de iteración
  error = 1
  iter = 0
  
  while ((error > tol) & (iter < maxiter)) {
    # E-step
    iter = iter + 1
    abilityold=ability
    parold=par
    z = cbind(matrix(1, q, 1), X) %*% t(par)
    plj = 1/(1 + exp(-1 * z))
    plj=cbind(plj,1-plj)
    L = matrix(1, s, q)
    # for (l in 1:s) for (k in 1:q) for (j in 1:n) L[l, k] = L[l, k] * (plj[k, j]^x[l, j]) * ((1 - plj[k, j])^(1 - x[l, j]))
    
    for (l in 1:s) 
      for (k in 1:q)
        L[l, k] = prod((plj[k,]^xx[l,]))
    # ----- Hay que cambiar esta linea. Va muy lenta ....
    
    Pl = L %*% A
    ability = matrix(0, s, dim)
    for (l in 1:s) for (j in 1:dim) {
      for (k in 1:q) ability[l, j] = ability[l, j] + X[k, j] * L[l, k] * A[k]
      ability[l, j] = ability[l, j]/Pl[l]
    }
    # M-step  -  Calculo de los parametros
    for (j in 1:n) par[j, ] = RidgeBinaryLogisticFit(x[, j], cbind(matrix(1, s, 1), ability), freq, tol, maxiter, penalization)
    
    lin = cbind(matrix(1, s, 1), ability) %*% t(par)
    PI = 1/(1 + exp(-1 * lin))
    logLik = sum(diagonal(freq) %*% log(PI^x * ((1 - PI)^(1 - x))))
    error = (logLik - logLikold)/abs(logLik)
    
    if (error < 0){
      ability=abilityold
      par=parold}
    else {
      print(c(iter, logLik, error))}
    logLikold = logLik
  }
  
  p=n
  Res=list()
  Res$Biplot="Binary Logistic"
  Res$RowCoordinates=ability
  Res$ColumnParameters=matrix(0,p,dimens+1)
  Res$NullDeviances=matrix(0,p,1)
  Res$ModelDeviances=matrix(0,p,1)
  Res$Deviances=matrix(0,p,1)
  Res$Dfs=matrix(0,p,1)
  Res$pvalues=matrix(0,p,1)
  Res$Bonferroni=matrix(0,p,1)
  Res$Nagelkerke=matrix(0,p,1)
  Res$R2=matrix(0,p,1)
  Res$PercentsCorrec=matrix(0,p,1)
  Res$DevianceTotal=0
  Res$p=1
  Res$TotalPercent=0
  Res$SSRes=matrix(0,p,1)
  Res$SSTot=matrix(0,p,1)
  
  for (i in 1:p){
    cat(paste(" ",i))
    y=x[,i]
    fit=RidgeBinaryLogistic(y,Res$RowCoordinates,tolerance = tol, maxiter = maxiter, penalization=penalization, cte=TRUE)
    Res$ColumnParameters[i,]=t(fit$beta)
    Res$ModelDeviances[i]=fit$Deviance
    Res$NullDeviances[i]=fit$NullDeviance
    Res$Deviances[i]=fit$Dif
    Res$Dfs[i]=fit$df
    Res$pvalues[i]=fit$p
    Res$Bonferroni[i]=(fit$p * p)* ((fit$p * p)<=1) + (((fit$p * p)>1))
    Res$R2[i]=fit$R2
    Res$CoxSnell[i]=fit$CoxSnell
    Res$Nagelkerke[i]=fit$Nagelkerke
    Res$MacFaden[i]=fit$MacFaden
    Res$PercentsCorrec[i]=fit$PercentCorrect
    Res$TotalPercent=Res$TotalPercent+sum(y==fit$Prediction)
    Res$SSRes[i]=fit$SSRes
    Res$SSTot[i]=fit$SSTot
  }
  
  d = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(dimens + 1)])^2))
  Res$Loadings = solve(diag(d)) %*% Res$ColumnParameters[, 2:(dimens + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)
  
  rownames(Res$ColumnParameters)=colnames(x)
  colnames(Res$ColumnParameters)=c("Const.",paste("Dim",1:dimens, sep=""))
  
  rownames(Res$ColumnParameters)=colnames(x)
  
  Res$TotalPercent=Res$TotalPercent/(s*p)
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$DevianceTotal=sum(Res$Deviances)
  
  Res$TotalSSRes=sum(Res$SSRes)
  Res$TotalSSTot=sum(Res$SSTot)
  
  nn=length(x)
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  Res$TotR2=1-(Res$TotalSSRes=sum(Res$SSRes)/Res$TotalSSTot)
  Res$TotalDf=sum(Res$Dfs)
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)
  Res$ClusterType="us"
  Res$Clusters = as.factor(matrix(1,nrow(Res$RowCoordinates), 1))
  Res$ClusterColors="blue"
  Res$ClusterNames="ClusterTotal"
  class(Res) = "Binary.Logistic.Biplot"
  
  return(Res)
  
}
