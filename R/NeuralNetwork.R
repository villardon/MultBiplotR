NeuralNetCost <- function(X, Y, hidden, theta){
  n=nrow(X)
  p=ncol(X)
  q=ncol(Y)
  sl=c(p, hidden, q)
  L=length(sl)
  Theta=list()
  ini=0
  for (i in 1 : (L-1)){
    elem=sl[i+1]*sl[i]
    end=ini+elem
    Theta[[i]]=matrix(theta[ini:end], nrow=sl[i+1])
    ini=end
  }
  
}


NeuralNetwork <- function(X, y, hidden=c(5,2,3), InitPar=NULL, rang=0.5, maxiter=100){

}