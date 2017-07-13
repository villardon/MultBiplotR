PCA.Bootstrap <- function(X, dimens=2, Scaling = "Standardize columns", B=1000, type="np"){
  # X is supposed to be transformed previously
  # If data are not standardized only non parametric Bootstrap should be used
  # For the moment, only reflection of the eigenvectors is considered
  # Some kind of rotation (Procrustes) may be needed 
  # type = c("np", "pa", "spper", "spres")
  # np: Non Parametric
  # pa: parametric (resuduals are suppose to have a normal distribution)
  # spper: Semi-parametric Residuals are permutated
  # spres: Semi-parametric Residuals are resampled
  
  if (is.data.frame(X)) X=as.matrix(X)
  n= dim(X)[1]
  p= dim(X)[2]
  
  rnames=rownames(X)
  cnames=colnames(X)
  dimnames=paste("Dim", 1:dimens)
  
  res=list()
  res$Type=type
  res$InitTransform=Scaling
  res$InitialData=X  
  XT=InitialTransform(X, transform = Scaling)$X
  res$TransformedDataData=XT
  res=list()
  res$Type=type
  res$X=X
  # Initial Calculation of the Principal Components
  acp=svd(XT, nu=dimens, nv=dimens)
  rownames(acp$u)=rnames
  colnames(acp$u)=dimnames
  rownames(acp$v)=cnames
  colnames(acp$v)=dimnames
  V=acp$v
  res$InitialSVD=acp
  # Initialization of the matrices of bootstrap estimates
  res$Samples=matrix(0, B, n)
  res$EigVal=matrix(0, B, min(n,p))
  res$Inertia=matrix(0, B, min(n,p))
  res$Us=array(0, dim=c(n, dimens, B))
  res$Vs=array(0, dim=c(p, dimens, B))
  res$As=array(0, dim=c(n, dimens, B))
  res$Bs=array(0, dim=c(p, dimens, B))
  res$Struct=array(0, dim=c(p, dimens, B))
  
  if (type=="np"){
    # Resampling individuals
    for (i in 1:B){
      samp=sample.int(n, size=n, replace=TRUE)
      res$Samples[i,]=samp
      XB=X[samp,]
      XB=InitialTransform(XB, transform = Scaling)$X
      acpB=svd(XB, nu=dimens, nv=dimens)
      signs=sign(diag(t(V) %*% acpB$v))
      acpB$v= acpB$v * matrix(1, p,1 ) %*% matrix(signs, ncol=dimens)
      acpB$u= acpB$u * matrix(1, n,1 ) %*% matrix(signs, ncol=dimens)
      rownames(acpB$u)=rnames
      colnames(acpB$u)=dimnames
      rownames(acpB$v)=cnames
      colnames(acpB$v)=dimnames
      
      res$EigVal[i,]=acpB$d
      res$Inertia[i,]=100*(acpB$d^2)/sum(acpB$d^2)
      res$Us[,,i]=acpB$u
      res$Vs[,,i]=acpB$v
      res$As[,,i]= XB %*% acpB$v
      res$Bs[,,i]= t(XB) %*% acpB$u
      res$Struct[,,i]=cor(XB, res$As[,,i])
    }
  }
  class(res)="PCA.Bootstrap"
  return(res)
}

plot.PCA.Bootstrap <- function(pcaboot, Eigenvalues=TRUE, Inertia=TRUE, EigenVectors=TRUE, Structure=TRUE, Squared=T){
  if (Eigenvalues){
    boxplot(pcaboot$EigVal^2, main="Bootstrap plot: Eigenvalues")}
  
  dimens=dim(pcaboot$Vs)[2]
  
  if (Inertia){
    boxplot(pcaboot$Inertia, main="Bootstrap plot: Percent of Explained Variance")}
  
  if (EigenVectors){
    op <- par(mfrow=c(dimens,1))
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Vs[,i,]))
      colnames(bb)
      boxplot(bb, main=paste("Bootstrap plot: Eigenvectors -", "PC",i))
      abline(0,0)
      }
    par(op)
  }
  
  
  if (Structure){
    op <- par(mfrow=c(dimens,1))
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]))
      colnames(bb)
      boxplot(bb, main=paste("Bootstrap plot: Correlations with the component -", "PC",i))
      abline(0,0)
    }
    par(op)
  }
  
  if (Squared){
    op <- par(mfrow=c(dimens,1))
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]^2))
      colnames(bb)
      boxplot(bb, main=paste("Bootstrap plot: Squared Correlations with the component -", "PC",i))
      abline(0,0)
    }
    par(op)
  }
  
}


