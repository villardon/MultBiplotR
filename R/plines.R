plines <- function(SetA,SetB,lin="dotted") {
  np<-nrow(SetA)
  for(i in 1:np) lines(rbind(SetA[i,1:2],SetB[i,1:2]),lty=lin) 
  return(NULL)
}



ones <- function(n) {
  if (length(n) == 1) 
    res = matrix(1, n, n)
  else res = array(1, n)
  return(res)
}



tr <- function(X) {
  tra = sum(diag(X))
  return(tra)
}