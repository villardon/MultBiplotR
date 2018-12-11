
sigmoide<-function(z){
  (1/(1+exp(-1*z)))
}

JLogBiplotReg <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n*r)+p*(r+1))], p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B[,-1]^2, na.rm = TRUE)/2
  return(J)
}


grLogBiplotReg <- function(par, X, r, lambda) { ## Gradient of 'fr'
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n*r)+p*(r+1))], p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradA=E%*%B[,-1]+lambda*A
  gradB=t(E)%*%cbind(rep(1,n),A)+lambda*cbind(rep(0,p),B[,-1])
  grad=c(c(gradA),c(gradB))
  return(grad)
}


# Loss Functions with LASSO penalization
JLogBiplotLASSO <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n*r)+p*(r+1))], p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(abs(par))
  return(J)
}


grLogBiplotLASSO <- function(par, X, r, lambda) { ## Gradient of 'fr'
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n*r)+p*(r+1))], p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradA=E%*%B[,-1]+lambda*sum(sign(A))
  gradB=t(E)%*%cbind(rep(1,n),A) + lambda*sign(cbind(rep(0,p),B[,-1]))
  grad=c(c(gradA),c(gradB))
  return(grad)
}


# Cost and gradients for the alternate algoritms
JLogBiplotRegB <- function(par, X, A, lambda) { # Cost to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]
  B=matrix(par, p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(B[,-1]^2, na.rm = TRUE)/2
  return(J)
}

JLogBiplotRegA <- function(par, X, B, lambda) { # Cost to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]-1
  A=matrix(par, n, r)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(A^2, na.rm = TRUE)/2
  return(J)
}

grLogBiplotRegB <- function(par, X, A, lambda) { ## Gradient to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]
  B=matrix(par, p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradB=t(E)%*%cbind(rep(1,n),A)+lambda*cbind(rep(0,p),B[,-1])
  grad=c(c(gradB))
  return(grad)
}


grLogBiplotRegA <- function(par, X, B, lambda) { ## Gradient to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]-1
  A=matrix(par, n, r)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradA=E%*%B[,-1]+lambda*A
  grad=c(c(gradA))
  return(grad)
}

