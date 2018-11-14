NIPALS.Biplot <- function(X, alpha = 1, dimension = 3, Scaling = 5, Type="Sparse", sup.rows = NULL, sup.cols = NULL, grouping=NULL, ...) {
  
  # Vamos a probar si esta cosa se actualiza
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

  if (Type=="Sparse") Biplot$Title = "Sparse NIPALS Biplot"
  if (Type=="Truncated") Biplot$Title = "Truncated NIPALS Biplot"
  if (Type=="Regular") Biplot$Title = "NIPALS Biplot"
    
  Biplot$Type = "PCA" 
  Biplot$call <- match.call()
  Biplot$Non_Scaled_Data = X
  Biplot$alpha=alpha
  Biplot$Dimension=dimension
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
  if (Scaling=="Within groups standardization") Biplot$Deviations = Data$ColStdDevs
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
  if (Type =="Sparse") SD = Sparse.NIPALSPCA(X, dimens = dimension, ...)
  if (Type =="Truncated") SD = Truncated.NIPALSPCA(X, dimens = dimension, ...)
  if (Type =="Regular") SD = NIPALSPCA(X, dimens = dimension, ...)
  
  
  a = SD$u %*% diag(SD$d)
  b = SD$v
  rownames(a) <- RowNames
  colnames(a) <- DimNames
  rownames(b) <- VarNames
  colnames(b) <- DimNames
  CorrXCP = cor(X, a, method = "pearson")
  # Relative contributions of the rows and columns
  sf = apply((X^2), 1, sum)
  cf=matrix(0,n,dimension)
  sc = apply((X^2), 2, sum)
  cfacum=matrix(0, n, dimension)
  ccacum=matrix(0, p, dimension)
  cf=matrix(0, n, dimension)
  cc=matrix(0, p, dimension)
  EV = matrix(0, dimension, 1)
  EVacum = matrix(0, dimension, 1)
  
  for (k in 1:dimension){
    Fitted.Val= a[,1:k] %*% t(b[,1:k])
    EVacum[k]=sum(Fitted.Val^2)/sum(X^2)
    cfacum[,k] = apply((Fitted.Val^2), 1, sum)/sf
    ccacum[,k] = apply((Fitted.Val^2), 2, sum)/sc
    if (k==1){
      EV[k]=EVacum[k]
      cc[,k]=ccacum[,k]
      cf[,k]=cfacum[,k]
    }
    else{
      EV[k]=EVacum[k] - EVacum[k-1]
      cc[,k]=ccacum[,k] - ccacum[,(k-1)]
      cf[,k]=cfacum[,k] - cfacum[,(k-1)]
    }
  }
  
  rownames(cfacum) = RowNames
  colnames(cfacum) = DimNames
  rownames(ccacum) = VarNames
  colnames(ccacum) = DimNames

  Inertia = round((EV/sum(X^2)) * 100, digits = 3)
  CumInertia = round((EVacum/sum(X^2)) * 100, digits = 3)
  
  if (alpha <= 1) {
    sca = sum(a^2)
    scb = sum(b^2)
    sca = sca/n
    scb = scb/p
    scf = sqrt(sqrt(scb/sca))
    a = a * scf
    b = b/scf
    if (nfs > 0) 
      as = as * scf
    if (ncs > 0) 
      bs = bs/scf
  } else {
    b = b %*% diag(SD$d)
    scf = 1
  }
  
  Biplot$nrows = n
  Biplot$ncols = p
  Biplot$nrowsSup = nfs
  Biplot$ncolsSup = ncs
  Biplot$dim = dimension
  Biplot$EigenValues = EV
  Biplot$Inertia = Inertia
  
  Biplot$CumInertia = CumInertia
  Biplot$V = SD$v
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