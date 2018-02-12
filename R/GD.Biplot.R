# Classical Biplot for Principal Components Analysis using Gradient Descent
GD.Biplot <- function(X, dimension = 2, Scaling = 5, lambda=0.01, OptimMethod="CG", Orthogonalize=FALSE, Algorithm="Alternated", sup.rows = NULL, 
                      sup.cols = NULL, grouping=NULL, tolerance=0.0001, num_max_iters=300, Initial="random") {
  if (is.data.frame(X)) 
    X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  # Setting the properties of data
  if (is.null(rownames(X))) 
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(X)
  if (is.null(colnames(X))) 
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  VarNames = colnames(X)
  
  DimNames = "Dim 1"
  for (i in 2:dimension) DimNames = c(DimNames, paste("Dim", i))
  if (is.null(sup.rows)) 
    nfs = 0
  else {
    if (!(p == ncol(sup.rows))) 
      stop("The #cols of the supplementary rows must be the same as the #cols of X")
    nfs = nrow(sup.rows)
    if (is.null(rownames(sup.rows))) 
      rownames(sup.rows) <- rownames(sup.rows, do.NULL = FALSE, prefix = "IS")
    colnames(sup.rows) <- VarNames
  }
  if (is.null(sup.cols)) 
    ncs = 0
  else {
    if (!(n == nrow(sup.cols))) 
      stop("The #rows of the supplementary columns must be the same as the #rows of X")
    ncs = ncol(sup.cols)
    if (is.null(colnames(sup.cols))) 
      colnames(sup.cols) <- colnames(sup.cols, do.NULL = FALSE, prefix = "SV")
    rownames(sup.cols) <- RowNames
  }
  Biplot = list()
  Biplot$Title = "Ridge Penalized Gradient Descent Biplot"
  Biplot$Type = "Alternated" 
  Biplot$call <- match.call()
  Biplot$Non_Scaled_Data = X
  Biplot$alpha="NA"
  Biplot$Dimension=dimension
  Biplot$Means = apply(X, 2, mean, na.rm=TRUE)
  Biplot$Medians = apply(X, 2, median, na.rm=TRUE)
  Biplot$Deviations = apply(X, 2, sd, na.rm=TRUE)
  Biplot$Minima = apply(X, 2, min, na.rm=TRUE)
  Biplot$Maxima = apply(X, 2, max, na.rm=TRUE)
  Biplot$P25 = apply(X, 2, quantile, na.rm=TRUE)[2, ]
  Biplot$P75 = apply(X, 2, quantile, na.rm=TRUE)[4, ]
  Biplot$GMean = mean(X, na.rm=TRUE)
  Biplot$Sup.Rows = sup.rows
  Biplot$Sup.Cols = sup.cols
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(Scaling)) 
    Scaling = ContinuousDataTransform[Scaling]
  Biplot$Initial_Transformation = Scaling
  Data = InitialTransform(X, sup.rows, sup.cols, transform = Scaling, grouping=grouping)
  X = Data$X
  rownames(X) = RowNames
  colnames(X) = VarNames
  Biplot$Scaled_Data = X
  Biplot$Scaled_Sup.Rows = Data$sup.rows
  Biplot$Scaled_Sup.Cols = Data$sup.cols
  if (nfs > 0) {
    rownames(Biplot$Scaled_Sup.Rows) <- rownames(sup.rows)
    colnames(Biplot$Scaled_Sup.Rows) <- colnames(sup.rows)
  }
  if (ncs > 0) {
    rownames(Biplot$Scaled_Sup.Cols) <- rownames(sup.cols)
    colnames(Biplot$Scaled_Sup.Cols) <- colnames(sup.cols)
  }
  # Calculating the Biplot
  
  
  # Initial (random) values for A and B
  r=dimension
  
  if (Initial=="random"){
    A=matrix(runif(n*r, -1, 1),n,r)
    B=matrix(runif(p*r, -1, 1), p, r)
  }
  
  if (Initial=="pca"){
    dvs=svd(X)
    A=dvs$u[,1:r] %*% diag(dvs$d[1:r])
    B=dvs$v[,1:r]
  }
  
  r=dimension
  
  if (Algorithm == "Alternated"){
    parA=c(c(A))
    parB=c(c(B))
    resbipB <- optimr(parB, fn=JBiplotRegB, gr=grBiplotRegB, method=OptimMethod, X=X, A=A, lambda=lambda)
    parB=resbipB$par
    J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
    err=1
    iter=0
    
    while( (err > tolerance) & (iter<num_max_iters)){
      iter=iter+1
      Jold=J
      #Update A
      resbipA <- optimr(parA, fn=JBiplotRegA, gr=grBiplotRegA, method=OptimMethod, X=X, B=B, lambda=lambda)
      parA=resbipA$par
      A=matrix(parA,n,r)
      if (Orthogonalize) {
        A=InitialTransform(A)$X
        A=Orthog(A)}
      #Update B
      resbipB <- optimr(parB, fn=JBiplotRegB, gr=grBiplotRegB, method=OptimMethod, X=X, A=A, lambda=lambda)
      parB=resbipB$par
      B=matrix(parB, p, r)
      J = sum((X - A %*% t(B))^2, na.rm = TRUE)/2 + lambda*sum(A^2, na.rm = TRUE)/2 + lambda*sum(B^2, na.rm = TRUE)/2
      err=(Jold-J)/Jold
      cat("\n",round(iter), round(J, 3), round(err,6))
    }
  }
  
  if (Algorithm == "Joint"){
    initpar=c(c(A),c(B))
    resbip <- optimr(initpar, fn=JBiplotReg, gr=grBiplotReg, method=OptimMethod, X=X, r=r, lambda=lambda)
    par=resbip$par
    A=matrix(par[1:(n*r)],n,r)
    B=matrix(par[(n*r+1):((n+p)*r)], p, r)
    if (Orthogonalize) {
      # A=InitialTransform(A)$X
      A=Orthog(A)
      parB=c(c(B))
      resbipB <- optimr(parB, fn=JBiplotRegB, gr=grBiplotRegB, method=OptimMethod, X=X, A=A, lambda=lambda)
      parB=resbipB$par
      B=matrix(parB, p, r)}
    
  }
  
  if (Algorithm == "Recursive"){
    XR=X
    for (i in 1:r){
      initpar=c(c(A[,i]),c(B[,i]))
      resbip <- optimr(initpar, fn=JBiplotReg, gr=grBiplotReg, method=OptimMethod, X=XR, r=1, lambda=lambda)
      par=resbip$par
      A1=matrix(par[1:n],n,1)
      B1=matrix(par[(n+1):(n+p)], p, 1)
      A[,i]=A1
      B[,i]=B1
      XR=XR- A1 %*% t(B1)
    }
    
  }
  
  
  # Matching the scales of row and column coordinates
  sca = sum(A^2)
  scb = sum(B^2)
  sca = sca/n
  scb = scb/p
  scf = sqrt(sqrt(scb/sca))
  A = A * scf
  B = B/scf
  
  
  scA=apply(A^2,2,sum)
  #A=A%*%diag(1/sqrt(scA))
  scB=apply(B^2,2,sum)
  #B=B%*%diag(1/sqrt(scB))
  EV=sqrt(scA)*sqrt(scB)
  Xh=A %*% t(B)
  X[which(is.na(X))]=Xh[which(is.na(X))]
  
  
  
  # Xh=A%*% d %*% t(B)
  sct= sum(X^2, na.rm = TRUE)
  CumInertia=matrix(0,nrow=r)
  cfacum=matrix(0,n,dimension)
  ccacum=matrix(0,p,dimension)
  for (i in 1:r){
    Xh=matrix(A[,1:i],nrow=n) %*% t(matrix(B[,1:i],nrow=p))
    cfacum[,i]= apply((Xh*X), 1, sum)^2 / (apply(X^2, 1, sum) * apply(Xh^2, 1, sum))
    ccacum[,i]= apply((Xh*X), 2, sum)^2 / (apply(X^2, 2, sum) * apply(Xh^2, 2, sum))
    E=X-Xh
    CumInertia[i]= 1-(sum(E^2)/sct)
  }
  
  ccacum=ccacum*100
  cfacum=cfacum*100
  cf = cfacum -cbind(rep(0,n),cfacum[,-r])
  cc = ccacum -cbind(rep(0,p),ccacum[,-r])
  CumInertia=CumInertia*100
  Inertia = CumInertia -c(0,CumInertia[-r])
  
  rownames(A) <- RowNames
  colnames(A) <- DimNames
  rownames(B) <- VarNames
  colnames(B) <- DimNames
  CorrXCP = cor(X, A, method = "pearson")
  
  rownames(cf) = RowNames
  colnames(cf) = DimNames
  rownames(cfacum) = RowNames
  colnames(cfacum) = DimNames
  # Relative contributions of the columns 
  rownames(cc) = VarNames
  colnames(cc) = DimNames
  rownames(ccacum) = VarNames
  colnames(ccacum) = DimNames
  
  Biplot$nrows = n
  Biplot$ncols = p
  Biplot$nrowsSup = nfs
  Biplot$ncolsSup = ncs
  Biplot$dim = dimension
  Biplot$EigenValues = EV
  Biplot$Inertia = Inertia
  
  Biplot$CumInertia = CumInertia
  Biplot$Structure = CorrXCP
  
  Biplot$RowCoordinates <- A
  Biplot$ColCoordinates <- B
  
  # Contributions
  Biplot$RowContributions <- cf
  Biplot$ColContributions <- cc
  Biplot$Scale_Factor = scf
  Biplot$ClusterType="us"
  Biplot$Clusters = as.factor(matrix(1,nrow(Biplot$RowContributions), 1))
  Biplot$ClusterColors="blue"
  Biplot$ClusterNames="Cluster 1"
  class(Biplot) <- "ContinuousBiplot"
  return(Biplot)
}

