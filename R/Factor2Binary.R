Factor2Binary <- function(y, Name=NULL){  
  if (is.null(Name)) Name="C-"
  ncat=length(levels(y))
  n=length(y)
  Z=matrix(0,n,ncat)
  for (i in 1:n)
    Z[i,as.numeric(y[i])]=1
  colnames(Z) <- paste(Name,levels(y),sep="")
  return(Z)
}