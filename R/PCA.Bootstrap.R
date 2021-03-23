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
  res$TransformedData=XT
  res$Type=type

  # Initial Calculation of the Principal Components
  acp=svd(XT, nu=dimens, nv=dimens)
  rownames(acp$u)=rnames
  colnames(acp$u)=dimnames
  rownames(acp$v)=cnames
  colnames(acp$v)=dimnames
  V=acp$v
  res$InitialSVD=acp
  res$InitScores=XT %*% acp$v
  res$InitCorr= cor(XT, res$InitScores)
  # Initialization of the matrices of bootstrap estimates
  res$Samples=matrix(0, B, n)
  res$EigVal=matrix(0, B, min(n,p))
  res$Inertia=matrix(0, B, min(n,p))
  res$Us=array(0, dim=c(n, dimens, B))
  res$Vs=array(0, dim=c(p, dimens, B))
  res$As=array(0, dim=c(n, dimens, B))
  res$Scores=array(0, dim=c(n, dimens, B))
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
      res$Scores[,,i]= XT %*% acpB$v
    }
  }
  
  if (type=="spres"){
    # Resampling residuals
    for (i in 1:B){
    }
  }
  
  #Confidence Intervals for the EigenValues
  res$MedEigVal=apply(res$EigVal^2, 2,mean)
  res$ICPercEigVal=apply(res$EigVal^2, 2,quantile, c(0.025, 0.975))
  sdev=apply(res$EigVal^2, 2,sd)
  res$ICBasEigVal=rbind(res$MedEigVal-1.96*sdev,res$MedEigVal+1.96*sdev)
 
  class(res)="PCA.Bootstrap"
  return(res)
}




