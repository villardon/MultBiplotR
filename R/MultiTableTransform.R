MultiTableTransform <- function(X, InitTransform="Standardize columns", dual=FALSE){
  ng = length(X) #Number of groups
  if (!dual){
    for (i in 1:ng){
      X[[i]] = TransformIni(X[[i]],transform = InitTransform)}
    return(X)
  }
  
  if (dual){
    Varnames=colnames(X[[1]])
    if (InitTransform=="Standardize columns"){
      nr=rep(0, ng)
      sd=NULL
      for (i in 1:ng){
        nr[i]=dim(X[[i]])[1]
        sd=rbind(sd,apply(X[[i]], 2, sd))
        X[[i]] = TransformIni(X[[i]],transform = "Column centering")
      }
      
      sds=sqrt(apply((diag(nr-1) %*% sd^2),2,sum)/sum(nr-1))
      for (i in 1:ng){
        X[[i]] = X[[i]] %*% diag(1/sds)
        colnames(X[[i]])=Varnames
      }
    }
  }
  else
    for (i in 1:ng){
      X[[i]] = TransformIni(X[[i]],InitTransform)}
  return(X)
}
