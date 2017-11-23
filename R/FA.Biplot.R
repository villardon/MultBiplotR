#  Biplot for Factor Analysis
FA.Biplot <- function(X, dimension = 3, Extraction="PC", Rotation="varimax", 
                      InitComunal="A1", normalize=FALSE, Scores= "Regression",  
                      MethodArgs=NULL, sup.rows = NULL, sup.cols = NULL, ...) {
  InitComunals=c("A1", "HSC", "MC")
  if (is.numeric(InitComunal)) 
    InitComunal = InitComunals[InitComunal]
  Scoress=c("Regression", "Bartlett")
  if (is.numeric(Scores)) 
    Scores = Scoress[Scores]
  
  Scaling=5;
  alpha=0;
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
  
  DimNames = paste("Factor", 1:dimension, sep="_")
  
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
  Biplot$Title = "Factor Analysis Biplot"
  Biplot$Type = "FA" 
  Biplot$call <- match.call()
  Biplot$alpha=alpha
  Biplot$Dimension=dimension
  Biplot$Non_Scaled_Data = X
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
  Biplot$GMean = mean(X)
  Biplot$Sup.Rows = sup.rows
  Biplot$Sup.Cols = sup.cols
  Extractions=c("PC", "IPF", "ML")
  if (is.numeric(Extraction)) 
    Extraction = Extractions[Extraction]
  Biplot$Initial_Transformation = Scaling
  Data = InitialTransform(X, sup.rows, sup.cols, transform = Scaling)
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
  R = cor(X)
  Biplot$R = R
  
  if (Extraction=="PC"){
    SD = svd(R)
    b = SD$u[,1:dimension] %*% diag(sqrt(SD$d[1:dimension]))
    SCT=p
    EV = SD$d[1:dimension]
  }
  
  if (Extraction=="IPF"){
    SD = PrinAxesFA(R, dimsol=dimension, method=InitComunal, tol=0.0001,  MaxIter=50)
    b = SD$A
    SCT=p
    EV = SD$Eigenvalues
  }
  
  if (Extraction=="ML"){
    SD = factanal(X, factors = dimension, rotation ="none")
    b = matrix(SD$loadings, p , dimension)
    rownames(b)=rownames(R)
    colnames(b)=paste("Factor",1:dimension, sep="_")
    SCT=p
    EV = apply(b^2, 2, sum)
  }
 # Rotaciones
  Rotations=c("quartimax", "varimax", "equamax", "parsimax", "entropy", "bifactorT", "oblimin", "quartimin", "oblimax", "simplimax", "bifactorQ" )
 if (is.numeric(Rotation)) 
   Rotation = Rotations[Rotation]

 if (is.element(Rotation, c("quartimax", "varimax", "equamax", "parsimax", "entropy", "bifactorT")))
   rotated=GPForth(b, normalize=normalize, method=Rotation, methodArgs=MethodArgs, ...)
 if (Rotation=="eaquamax")
   rotated=cfT(b,kappa=dimension/(2*p))
 if (Rotation=="parsimax")
   rotated=cfT(b,kappa=(dimension-1)/(p+dimension-2))
 if (is.element(Rotation, c( "oblimin", "quartimin", "oblimax", "simplimax", "bifactorQ" )))
   rotated=GPFoblq(b, normalize=normalize, method=Rotation, methodArgs=MethodArgs, ...)
 if (Rotation != "none")
  b=rotated$loadings

  # Me falta calcular la rotación, hay un enorme lío y bastante diferencia con Matlab
  
  h=apply(b^2,1,sum)
  u=1-h
  
  switch(Scores, Regression = {
    a= X %*% diag(1/u^2) %*% b %*% solve(t(b) %*% diag(1/u^2) %*% b)
  }, Bartlett={
    a=X %*% solve(R) %*% b
  })
  
  Inertia = round((EV/SCT) * 100, digits = 3)
  CumInertia = cumsum(Inertia)
  Biplot$Communalities=apply(b^2,1,sum)
  Biplot$Uniqueness=1-Biplot$Communalities
  rownames(a) <- RowNames
  colnames(a) <- DimNames
  rownames(b) <- VarNames
  colnames(b) <- DimNames
  CorrXCP = cor(X, a, method = "pearson")
  # Relative contributions of the rows
  
  sf = apply((X^2), 1, sum)
  cf=matrix(0,n,dimension)
  for (k in 1:dimension)
    cf[,k]= round((a[,k]*sqrt(EV[k]))^2/ sf*100, digits = 2)
  rownames(cf) = RowNames
  colnames(cf) = DimNames
  cfacum = t(apply(cf, 1, cumsum))
  # Relative contributions of the rows 
  
  cc=round((b^2)*100, digits = 2)
  rownames(cc) = VarNames
  colnames(cc) = DimNames
  ccacum = t(apply(cc, 1, cumsum))
  if (nfs > 0) {
    as = sup.rows %*% SD$v
    sf = apply((sup.rows^2), 1, sum)
    cfs = round((solve(diag(sf)) %*% as^2) * 100, digits = 2)
    rownames(cfs) = rownames(sup.rows)
    colnames(cfs) = DimNames
  }
  if (ncs > 0) {
    bs = t(sup.cols) %*% SD$u
    sc = apply((sup.cols^2), 2, sum)
    ccs = round((solve(diag(sc)) %*% bs^2) * 100, digits = 2)
    rownames(ccs) = colnames(sup.cols)
    colnames(ccs) = DimNames
  }

  scf = 1
  
  Biplot$nrows = n
  Biplot$ncols = p
  Biplot$nrowsSup = nfs
  Biplot$ncolsSup = ncs
  Biplot$dim = dimension
  Biplot$EigenValues = EV
  Biplot$Inertia = Inertia
  
  Biplot$CumInertia = CumInertia

  Biplot$Structure = CorrXCP
  if (nfs == 0) 
    Biplot$RowCoordinates <- a
  else Biplot$RowCoordinates <- rbind(a, as)
  if (ncs == 0) 
    Biplot$ColCoordinates <- b
  else Biplot$ColCoordinates <- rbind(b, bs)
  
  # Contributions
  if (nfs == 0) 
    Biplot$RowContributions <- cf
  else Biplot$RowContributions <- rbind(cf, cfs)
  if (ncs == 0) 
    Biplot$ColContributions <- cc
  else Biplot$ColContributions <- rbind(cc, ccs)
  Biplot$Scale_Factor = scf
  Biplot$ClusterType="us"
  Biplot$Clusters = as.factor(matrix(1,nrow(Biplot$RowContributions), 1))
  Biplot$ClusterColors="blue"
  Biplot$ClusterNames="Cluster 1"
  class(Biplot) <- "ContinuousBiplot"
  return(Biplot)
}

