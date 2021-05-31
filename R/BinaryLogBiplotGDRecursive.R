# Logistic Biplot with Gradient Descent
BinaryLogBiplotGDRecursive <- function(X, freq = matrix(1, nrow(X), 1),  dim = 2, tolerance = 1e-04,  penalization=0.2,
                              num_max_iters=100, RotVarimax = FALSE, OptimMethod="CG", Initial="random", ...) {
  # Joint algorithm for logistic biplots
  X=as.matrix(X)
  indnames=rownames(X)
  varnames=colnames(X)
  n <- nrow(X)
  p <- ncol(X)
  # Estimation of parameters A and B
  r=dim
  # Fitting the constants
  b0=matrix(0,nrow=p, ncol=1)
  for (j in 1:p)
    b0[j]=RidgeBinaryLogistic(y=X[,j], matrix(1,n,1), penalization = penalization)$beta
  B=b0
  A=matrix(1,nrow=n, ncol=1)
  H=sigmoide(A %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + penalization*sum(B[,-1]^2, na.rm = TRUE)/2
  ## Fitting the rest of the dimensions recursively
  for (k in 1:r){
    parA=X[,k]
    parB=rep(1,p)/sqrt(p)
    A=cbind(A,parA)
    B=cbind(B,parB)
    H=sigmoide(A %*% t(B))
    J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + penalization*sum(B[,-1]^2, na.rm = TRUE)/2
    err=1
    iter=0
    while( (err > tolerance) & (iter<num_max_iters)){
      iter=iter+1
      Jold=J
      #Update B
      resbipB <- optimr(parB, fn=JLogBiplotRegBRec, gr=grLogBiplotRegBRec, method=OptimMethod, X=X, A=A, B=B,lambda=penalization)
      parB=resbipB$par
      #parB=parB/sqrt(sum(parB^2))
      B[,k+1]=parB
      #Update A
      resbipA <- optimr(parA, fn=JLogBiplotRegARec, gr=grLogBiplotRegARec, method=OptimMethod, X=X, A=A, B=B, lambda=penalization)
      parA=resbipA$par
      parA=(parA-mean(parA))/sd(parA)
      A[,k+1]=parA
      H=sigmoide(A %*% t(B))
      J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2 + penalization*sum(B[,-1]^2, na.rm = TRUE)/2
      err=abs(Jold-J)/J
      cat("\n",round(iter), round(J, 3), round(err,6))
    }
  }
  A=A[,-1]
  if (RotVarimax) {
    BB = varimax(B[, 2:(dim + 1)], normalize = FALSE)
    A = A %*% BB$rotmat
    B[, 2:(dim + 1)] = B[, 2:(dim + 1)] %*% BB$rotmat
  }
  rownames(A)=indnames
  colnames(A)=paste("dim",1:dim)
  rownames(B)=varnames
  colnames(B)=c("intercept", paste("dim",1:dim))
  Res=list()
  Res$Data=X
  Res$Dimension=dim
  Res$Penalization=penalization
  Res$Tolerance=tolerance
  Res$OptimMethod=OptimMethod
  Res$Biplot="Binary Logistic (Recursive Gradient Descent)"
  Res$Type= "Binary Logistic (Recursive Gradient Descent)"
  Res$RowCoordinates=A
  #Res$ColumnParameters=cbind(b0,B)
  Res$ColumnParameters=B
  esp = cbind(rep(1,n), A) %*% t(B)
  pred = exp(esp)/(1 + exp(esp))
  pred2 = matrix(as.numeric(pred > 0.5), n, p)
  acier = matrix(as.numeric(round(X) == pred2), n, p)
  acierfil = 100*apply(acier,1,sum)/p
  aciercol = 100*apply(acier,2,sum)/n
  presences=apply(X, 2, sum)
  absences=n-presences
  sens = apply((acier==1) & (X==1), 2, sum)/presences
  spec = apply((acier==1) & (X==0), 2, sum)/absences
  totsens = sum((acier==1) & (X==1))/sum(presences)
  totspec = sum((acier==1) & (X==0))/sum(absences)
  gfit = (sum(sum(acier))/(n * p)) * 100
  esp0 = matrix(rep(1,n), n,1) %*% B[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))
  d1 = -2 * apply(X * log(pred0) + (1 - X) * log(1 - pred0),2,sum)
  d2 = -2 * apply(X * log(pred) + (1 - X) * log(1 - pred),2,sum)
  d = d1 - d2
  ps = matrix(0, p, 1)
  for (j in 1:p) ps[j] = 1 - pchisq(d[j], 1)
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
  Res$Loadings = diag(1/dd) %*% Res$ColumnParameters[, 2:(dim + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)
  nn=n*p
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  Res$R2 = apply((X-H)^2,2, sum)/apply((X)^2,2, sum)
  Res$TotR2 = sum((X-H)^2) /sum((X)^2)
  pred= matrix(as.numeric(H>0.5),n , p)
  verdad = matrix(as.numeric(X==pred),n , p)
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
  return(Res)
}


