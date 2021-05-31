JLogBiplotRegRec <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n*r)+p*(r+1))], p, r+1)
  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B[,-1]^2, na.rm = TRUE)/2
  return(J)
}


# Cost and gradients for the alternate algoritms
JLogBiplotRegBRec <- function(par, X, A, B, lambda) { # Cost to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]-1
  B[, r+1]=par
  H=sigmoide(A%*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(B[,-1]^2, na.rm = TRUE)/2
  return(J)
}

JLogBiplotRegARec <- function(par, X, A, B, lambda) { # Cost to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]-1
  A[,r+1]=par
  H=sigmoide(A %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + lambda*sum(A^2, na.rm = TRUE)/2
  return(J)
}

grLogBiplotRegBRec <- function(par, X, A, B, lambda) { ## Gradient to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]-1
  B[, r+1]=par
  H=sigmoide(A %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradB=t(E)%*%A+lambda*cbind(rep(0,p),B[,-1])
  grad=c(c(gradB[,r+1]))
  return(grad)
}


grLogBiplotRegARec <- function(par, X, A, B, lambda) { ## Gradient to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]-1
  A[,r+1]=par
  H=sigmoide(A %*% t(B))
  E = H-X
  E[which(is.na(X))]=0
  gradA=E%*%B+lambda*A
  grad=c(c(gradA[,r+1]))
  return(grad)
}

