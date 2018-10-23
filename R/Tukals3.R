Tuckals3 <- function(X, Scaling = 5, A=NULL, B=NULL, C=NULL , P=2, Q=2, R=2, tolerance=0.000001, maxiter=100){
  # Por el momento la transformación inicial será para las variables que, generalmente, se colocan en el segundo modo
  # Transformaremos cada variable en cada ocasión, es decir, para cada combinación variable-ocasión centramos sobre los individuos
  # Dependiendo de la transformación los biplots obtenidos pueden ser distintos
  
  K = length(X)
  I = dim(X[[1]])[1]
  J = dim(X[[1]])[2]
  
  Inames=rownames(X[[1]])
  Jnames=
  for (k in 1:K){
    X[[k]]=InitialTransform(X[[k]], transform = Scaling)$X
  }
  
  X1 = X[[1]]
  for (i in 2:K)
    X1 = cbind(X1, X[[i]])
  
  X2 = t(X[[1]])
  for (j in 2:K)
    X2 = cbind(X2, t(X[[j]]))
  
  X3 = matrix(0, K, I*J)
  
  for (k in 1:K)
    X3[k,]=vec(X[[k]])$v
  
  # Initial values for A, B, C y G
  if (is.null(A)) A=svd(X1)$u[,1:P]
  if (is.null(B)) B=svd(X2)$u[,1:Q]
  if (is.null(C)) C=svd(X3)$u[,1:R]
  
  K1=kronecker(C,B)
  
  G=t(A) %*% X1 %*% K1
  
  Xest=A %*% G %*% t(K1)
  
  SCRold = sum((X1-Xest)^2)
  error=1
  iter=0
  while ((error>tolerance) & (iter<maxiter)){
    iter=iter+1
    A=svd(X1 %*% K1)$u[,1:P]
    K2=kronecker(C,A)
    B=svd(X2 %*% K2)$u[,1:Q]
    K3=kronecker(A, B)
    C=svd(X3 %*% K3)$u[,1:R]
    K1=kronecker(C,B)
    G=t(A) %*% X1 %*% K1
    Xest=A %*% G %*% t(K1)
    SCRnew = sum((X1-Xest)^2)
    error=abs(SCRold-SCRnew)
    SCRold=SCRnew
    cat("Iteration :", iter, " -  Error:",error, "\n")
  }
  
  SCE=sum(G^2)
  SCT=sum(X1^2)
  ExPer=100*SCE/SCT
  
  Result=list(X=X, A=A, B=B, C=C, G=G, RSS=SCRnew, ESS=SCE, TSS=SCT, PercentExplained=ExPer)
  class(Result)="Tuckals3"
  return(Result)
}



