CrissCross <- function(x,w=matrix(1,dim(x)[1],dim(x)[2]),dimens=2, a0=NULL, b0=NULL, maxiter=100, tol=10e-5, addsvd=TRUE, lambda=0){
  n=dim(x)[1]
  p=dim(x)[2]
  x[which(w==0)]=0
  DimNames = "Dim 1"
  for (i in 2:dimens) DimNames = c(DimNames, paste("Dim", i))
  
  if (is.null(a0))
    a0=matrix(rnorm(n*dimens),n,dimens)
  if (is.null(b0))
    b0=matrix(rnorm(p*dimens),p,dimens)
  
  esp=(a0 %*% t(b0))
  resid=(x-esp)
  gfit0= sum((resid*w)^2)/sum((x*w)^2)
  error=0.1
  iter=0
  a=a0
  b=b0
  while ((error>tol) & (iter<maxiter)){
    iter=iter+1
    for (i in 1:n)
      a[i,]= ginv(t(b0) %*% diag(w[i,]+lambda) %*% b0) %*% t(b0) %*% diag(w[i,]) %*% (x[i,])
    
    for (j in 1:p)
    b[j,]=ginv(t(a) %*% diag(w[,j]+lambda) %*% a) %*% t(a) %*% diag(w[,j]) %*% (x[,j])
    
    esp=(a %*% t(b))
    resid=(x-esp)
    gfit= sum((resid*w)^2)/sum((x*w)^2)
    error = abs(gfit0-gfit)
    a0=a
    b0=b
    gfit0=gfit
    print(c(iter,error))
  }
  
  if (addsvd){
    sol=svd(esp)
    a=sol$u[,1:dimens] %*% sqrt(diag(sol$d[1:dimens]))
    b=sol$v[,1:dimens] %*% sqrt(diag(sol$d[1:dimens]))
  }
  
  
  rowfit=matrix(0,n,dimens)
  colfit=matrix(0,p,dimens)
  dimfit=matrix(0,dimens,1)
  rownames(a)=rownames(x)
  colnames(a)=DimNames
  rownames(b)=colnames(x)
  colnames(b)=DimNames
  
  for (k in 1:dimens){
    esp=(a[,k] %*% t(b[,k]))
    for (i in 1:n)
    rowfit[i,k]= wcor(x[i,], esp[i,], w[i,])^2
  for (j in 1:p)
    colfit[j,k]= wcor(x[,j], esp[,j], w[,j])^2
  
  }
  
  rownames(rowfit)=rownames(x)
  rownames(colfit)=colnames(x)
  colnames(rowfit)=DimNames
  colnames(colfit)=DimNames
  
  gfit=1-gfit
  
  CrissCross=list()
  CrissCross$Title = "ALS Biplot (CrissCross)"
  CrissCross$Weigths=w
    CrissCross$nrows=n
    CrissCross$ncols=p
    CrissCross$nrowsSup=0
    CrissCross$ncolsSup=0
    CrissCross$Dimension=dimens
    CrissCross$EigenValues
    CrissCross$Inertia
    CrissCross$CumInertia
    CrissCross$Structure
  CrissCross$RowCoordinates=a
  CrissCross$ColCoordinates=b
  CrissCross$RowContributions=rowfit*100
  CrissCross$ColContributions=colfit*100
  CrissCross$Scale_Factor=1
  CrissCross$alpha=1
  CrissCross$Dimension=dimens
  class(CrissCross) <- "ContinuousBiplot"
  return(CrissCross)
  
}

