

vec  <- function(X){
  n=dim(X)[1]
  p=dim(X)[2]
  z=X[,1]
  for (i in 2:p)
    z=c(z,X[,i])
  z
}


JUnfolding <- function(par, D, r) {
  n=dim(D)[1]
  p=dim(D)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  DD = DistUnfold(A, B)
  J = sum((D - DD)^2, na.rm = TRUE)/2
  return(J)
}

grUnfolding <- function(par, D, r) { ## Gradient of 'fr'
  n=dim(D)[1]
  p=dim(D)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  DD =DistUnfold(A, B)
  gradA=matrix(0,n,r)
  gradB=matrix(0,p,r)
  for (i in 1:r){
    gradA[,i] = apply((matrix(A[,i],n,1) %*% matrix(1,1,p)- matrix(1,n,1) %*% matrix(B[,i],1,p))*(1-D/DD), 1, sum) * A[,i]
    gradB[,i] = apply((matrix(1,n,1) %*% matrix(B[,i],1,p) - matrix(A[,i],n,1) %*% matrix(1,1,p))*(1-D/DD), 2, sum) * B[,i]
  }
  grad=c(c(gradA),c(gradB))
  return(grad)
}

JUnfoldingReg <- function(par, D, r, lambda) {
  n=dim(D)[1]
  p=dim(D)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  DD = DistUnfold(A, B)
  J = sum((D^2 - DD^2)^2, na.rm = TRUE)/4 + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
  return(J)
}

grUnfoldingReg <- function(par, D, r, lambda) { ## Gradient of 'fr'
  n=dim(D)[1]
  p=dim(D)[2]
  A=matrix(par[1:n*r],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  E = D^2 - DistUnfold(A, B)^2
  gradA = sum(E)*A - E%*%B + lambda*A
  gradB = sum(E)*B - t(E)%*%A +lambda*B
  grad=c(c(gradA),c(gradB))
  return(grad)
}


JUnfoldingLASSO <- function(par, X, r, lambda) {
  n=dim(X)[1]
  p=dim(X)[2]
  A=matrix(par[1:(n*r)],n,r)
  B=matrix(par[(n*r+1):((n+p)*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(abs(A), na.rm = TRUE) + lambda*sum(abs(B), na.rm = TRUE)
  return(J)
}

grUnfoldingLASSO <- function(par, X, r, lambda) { ## Gradient of 'fr'
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
JUnfoldingRegB <- function(par, X, A, lambda) { # Cost to estimate B
  n=dim(X)[1]
  p=dim(X)[2]
  r=dim(A)[2]
  B=matrix(par[1:(p*r)], p, r)
  J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
  return(J)
}

JUnfoldingRegA <- function(par, X, B, lambda) { # Cost to estimate A
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

grUnfoldingRegA <- function(par, X, B, lambda) { ## Gradient to estimate A
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



# Unfolding based on gradient descent
UnfoldingGD <- function(P, W = matrix(1, dim(P)[1], dim(P)[2]), Constrained = FALSE, ENV = NULL, model = "Ratio", condition = "Matrix", OptimMethod="CG", r = 2, maxiter = 100, 
                            tolerance = 1e-05, lambda = 0.01, omega = 0, X0, Y0) {
  # Using Gradient Descent
  n = dim(P)[1]
  m = dim(P)[2]
  DimNames = "Dim-1"
  for (i in 2:r) DimNames = c(DimNames, paste("Dim-", i, sep = ""))
  
  ColNames = colnames(P)
  RowNames = rownames(P)
  
  # Calculo de las distancias en la configuracion
  D = DistUnfold(X0, Y0)
  
  # Calculo de las disparidades optimas
  # [E,pstress]=cvpenalize(D,P,W,lambda,omega);
  
  DH = dhatsunfol(P, D, W, modelo = model, condicion = condition)
  pstress = sum(sum((W * (D - DH$Dh))^2))
  errorest = 1
  k = 0
  history = NULL
  print(lambda)
  while ((k <= maxiter) & (errorest > tolerance)) {
    k = k + 1
    initpar=c(c(X0),c(Y0))
    update <- optimr(initpar, fn=JUnfolding, gr=grUnfolding, method=OptimMethod, D=DH$Dh, r=r) 
    par=update$par
    X=matrix(par[1:(n*r)],n,r)
    Y=matrix(par[(n*r+1):((n+m)*r)], m, r)
    D = DistUnfold(X, Y)
    #[E,pstress1]=cvpenalize(D,Dh,W,lambda,omega)
    DH = dhatsunfol(P, D, W, model, condition)
    newpstress = sum(sum((W * (D - DH$Dh))^2))
    errorest = (pstress - newpstress)^2
    pstress = newpstress
    X0=X
    Y0=Y
    history = rbind(history, c(k, errorest))
    print(c(k, errorest))
  }
  
  rownames(X0) = RowNames
  colnames(X0) = DimNames
  rownames(Y0) = ColNames
  colnames(Y0) = DimNames
  
  dmean = sum(sum(D * W))/sum(sum(W))
  stress1 = sum(sum(((D - DH$Dh)^2) * W))/sum(sum(((D)^2) * W))
  stress2 = sum(sum(((D - DH$Dh)^2) * W))/sum(sum(((D - dmean)^2) * W))
  rs = cor(as.vector(DH$Dh * W), as.vector(D * W))
  rsq = rs^2
  
  sstress1 = sqrt(sum(sum(((D^2 - DH$Dh^2)^2) * W))/sum(sum(((D^2)^2) * W)))
  dmean2 = sum(sum((D^2) * W))/sum(sum(W))
  sstress2 = sqrt(sum(sum(((D^2 - DH$Dh^2)^2) * W))/sum(sum(((D^2 - dmean2)^2) * W)))
  rho = cor(as.vector(DH$Dh * W), as.vector(D * W), method = "spearman")
  tau = cor(as.vector(DH$Dh * W), as.vector(D * W), method = "kendall")
  
  fitcols = matrix(0, m, 1)
  for (i in 1:m) fitcols[i] = cor((DH$Dh * W)[, i], (D * W)[, i])^2
  rownames(fitcols) = ColNames
  colnames(fitcols) = "R-Squared"
  
  colnames(history) = c("Iteration", "Error")
  rownames(DH$Tol) = colnames(P)
  colnames(DH$Tol) = "Tolerance"
  if (Constrained) 
    Analysis = "Constrained Unfolding"
  else Analysis = "Unfolding"
  
  Unfold = list()
  Unfold$call <- match.call()
  Unfold$Analysis = Analysis
  Unfold$nsites = n
  Unfold$nspecies = m
  Unfold$nvars = NULL
  Unfold$alpha = NULL
  Unfold$dimens = r
  Unfold$Abundances = P
  Unfold$TransformedAbundances = P
  Unfold$Disparities = DH$Dh
  Unfold$Distances = D
  Unfold$Environment = ENV
  Unfold$TransfEnvironment = ENV
  Unfold$Weighted_Averages = NULL
  Unfold$Minima_Z = NULL
  Unfold$Maxima_Z = NULL
  Unfold$Medians_Z = NULL
  Unfold$P25_Z = NULL
  Unfold$P75_Z = NULL
  Unfold$Means = NULL
  Unfold$Deviations = NULL
  Unfold$IterationHistory = history
  Unfold$X = X0
  Unfold$Y = Y0
  Unfold$Tolerance = DH$Tol
  Unfold$RawStress = pstress
  Unfold$Stress1 = stress1
  Unfold$Stress2 = stress2
  Unfold$Sstress1 = sstress1
  Unfold$Sstress2 = sstress2
  Unfold$Pearson = rs
  Unfold$Spearman = rho
  Unfold$Kendall = tau
  Unfold$RSQ = rsq
  Unfold$ColumnsFit = fitcols
  Unfold$model = model
  Unfold$condition = condition
  class(Unfold) = "Unfolding"
  return(Unfold)
}