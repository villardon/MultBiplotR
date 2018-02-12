JBiplot <- function(par, X, r) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2
  return(J)
}

grBiplot <- function(par, X, r) { ## Gradient of 'fr'
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:n*r],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  E = A %*% t(B) - X
  E[which(is.na(X))]=0
  gradA=E%*%B
  gradB=t(E)%*%A
  grad=c(c(gradA),c(gradB))
  return(grad)
}


JBiplotReg <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
  return(J)
}


grBiplotReg <- function(par, X, r, lambda) { ## Gradient of 'fr'
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  E = A %*% t(B) - X
  E[which(is.na(X))]=0
  gradA=E%*%B+lambda*A
  gradB=t(E)%*%A+lambda*B
  grad=c(c(gradA),c(gradB))
  return(grad)
}



JBiplotLASSO <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(abs(A), na.rm = TRUE) + lambda*sum(abs(B), na.rm = TRUE)
  return(J)
}

grBiplotLASSO <- function(par, X, r, lambda) { ## Gradient of 'fr'
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  E = A %*% t(B) - X
  E[which(is.na(X))]=0
  gradA=E%*%B+lambda*sum(sign(A))
  gradB=t(E)%*%A+lambda* sum(sign(B))
  grad=c(c(gradA),c(gradB))
  return(grad)
}



# Cost and gradients for the alternate algoritms
JBiplotRegB <- function(par, X, A, lambda) { # Cost to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]
  B=matrix(par[1:(p*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
  return(J)
}

JBiplotRegA <- function(par, X, B, lambda) { # Cost to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]
  A=matrix(par[1:(n*r)],n,r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(A^2, na.rm = TRUE)/2 
  return(J)
}

grBiplotRegB <- function(par, X, A, lambda) { ## Gradient to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]
  B=matrix(par[1:(p*r)], p, r)
  E = A %*% t(B) - X
  E[which(is.na(X))]=0
  gradB=t(E)%*%A+lambda*B
  grad=c(c(gradB))
  return(grad)
}

grBiplotRegA <- function(par, X, B, lambda) { ## Gradient to estimate A
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(B)[2]
  A=matrix(par[1:(n*r)],n,r)
  E = A %*% t(B) - X
  E[which(is.na(X))]=0
  gradA=E%*%B+lambda*A
  grad=c(c(gradA))
  return(grad)
}

