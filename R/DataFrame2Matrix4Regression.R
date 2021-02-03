# The procedure keeps the numerical variables and convert factors to dummy variables for regression

DataFrame2Matrix4Regression <- function(X, last=TRUE, Intercept=FALSE){
  if (!is.data.frame(X)) stop("You must provide a data frame to prepare for regression")
  n=dim(X)[1]
  p=dim(X)[2]
  rnames=rownames(X)
  names=colnames(X)
  print(names)
  XN=NULL
  newnames=NULL
  for (i in 1:p){
    if (is.numeric(X[[i]])) {
      XN=cbind(XN,X[[i]])
      newnames=c(newnames, names[i])
    }
    if (is.factor(X[[i]])){
      Z=Factor2Binary(X[[i]], Name=names[i])
      pp=dim(Z)[2]
      nn=colnames(Z)
      if (last) {Z=as.matrix(Z[,-pp])
      nn=nn[-pp]}
      else {Z=as.matrix(Z[,-1])
            nn=nn[-1]}
      newnames=c(newnames, nn)
      XN=cbind(XN,Z)
    }
  }

  if (Intercept){
    XN=cbind(rep(1,n),XN)
    newnames=c("Intercept", newnames)
  }
  colnames(XN)<- newnames
  rownames(XN) =rnames
  return(XN)
}