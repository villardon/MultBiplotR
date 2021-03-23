AddIndiv2LogBiplotGD <- function(LB, X, freq = matrix(1, nrow(X), 1), tolerance = 1e-04,  penalization=0.01, 
                              num_max_iters=100, RotVarimax = FALSE, seed = 0, OptimMethod="CG", Initial="random",
                              Orthogonalize=FALSE, Algorithm = "Joint", ...) {
  
  
  X=as.matrix(X)
  indnames=rownames(X)
  varnames=colnames(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Estimation of parameters A and B
  r=LB$Dimension
  if (!is.null(LB$Seed)) set.seed(LB$Seed)
  
  B=LB$ColumnParameters
  A=X %*% B[,-1]
  parA=vec(A)$v
  
  
  resbipA <- optimr(parA, fn=JLogBiplotRegA, gr=grLogBiplotRegA, method=LB$OptimMethod, X=X, B=B, lambda=LB$Penalization)
  parA=resbipA$par
  A=matrix(parA,n,r)
  rownames(A)=rownames(X)

  H=sigmoide(cbind(rep(1,n),A) %*% t(B))
  J=sum(-1*X*log(H)-(1-X)*log(1-H), na.rm = TRUE)/2  + LB$Penalization*sum(A^2, na.rm = TRUE)/2 + LB$Penalization*sum(B[,-1]^2, na.rm = TRUE)/2
  
  
}
  