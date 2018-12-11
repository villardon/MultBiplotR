CostBinPLS1A <- function(par, Y, U, Q) { # Cost to estimate U
  if (is.vector(U)) U=matrix(U, ncol=1)
  if (is.vector(Q)) Q=matrix(Q, nrow=1)
  n=dim(Y)[1]
  p=dim(Y)[2]
  r=dim(Q)[2]
  H=sigmoide(U %*% t(Q))
  J=sum(-1*Y*log(H)-(1-Y)*log(1-H), na.rm = TRUE)/2
  return(J)
}

grBinPLS1A <- function(par, Y, U, Q) { ## Gradient to estimate U
  if (is.vector(U)) U=matrix(U, ncol=1)
  if (is.vector(Q)) Q=matrix(Q, nrow=1)
  n=dim(Y)[1]
  p=dim(Y)[2]
  r=dim(Q)[2]
  H=sigmoide(U %*% t(Q))
  E = H-Y
  E[which(is.na(Y))]=0
  gradA=E%*%Q[,r]
  grad=c(c(gradA))
  return(grad)
}

CostBinPLS1B <- function(par, Y, U, Q) { # Cost to estimate Q
  if (is.vector(U)) U=matrix(U, ncol=1)
  if (is.vector(Q)) Q=matrix(Q, nrow=1)
  n=dim(Y)[1]
  p=dim(Y)[2]
  r=dim(U)[2]
  H=sigmoide(U %*% t(Q))
  J=sum(-1*Y*log(H)-(1-Y)*log(1-H), na.rm = TRUE)/2
  return(J)
}

grBinPLS1B <- function(par, Y, U, Q) { ## Gradient to estimate B
  if (is.vector(U)) U=matrix(U, ncol=1)
  if (is.vector(Q)) Q=matrix(Q, nrow=1)
  n=dim(Y)[1]
  p=dim(Y)[2]
  r=dim(U)[2]
  H=sigmoide(U %*% t(Q))
  E = H-Y
  E[which(is.na(Y))]=0
  gradB=t(E)%*% U[,r]
  grad=c(c(gradB))
  return(grad)
}



