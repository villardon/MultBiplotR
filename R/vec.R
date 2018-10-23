vec <- function(X, ByColumns=TRUE){
  n=dim(X)[1]
  p=dim(X)[2]
  
  if (ByColumns){
    v=numeric()
    for (j in 1:p)
      v=c(v,X[,j])
  }
  return(list(v=v, ByColumns=ByColumns, n=n, p=p))
}