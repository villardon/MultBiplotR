# Classical Biplot for Principal Components Analysis
WPCA.Biplot <- function(X, CM=NULL, RM=NULL, alpha = 1, dimension = 3, Scaling = 5, sup.rows = NULL, sup.cols = NULL, grouping=NULL) {
  # Vamos a probar si esta cosa se actualiza
  if (is.data.frame(X)) 
    X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  
  if (is.null(CM)) CM=diag(n)
  if (is.null(RM)) RM=diag(p)

  # Setting the properties of data
  if (is.null(rownames(X))) 
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(X)
  if (is.null(colnames(X))) 
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  VarNames = colnames(X)
  
  rownames(CM)=rownames(X)
  colnames(CM)=rownames(X)
  rownames(RM)=colnames(X)
  colnames(RM)=colnames(X)
  
  DimNames = paste("Dim", 1:dimension)
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
  Biplot$Title = "GWPCA Biplot"
  Biplot$Type = "GWPCA" 
  Biplot$call <- match.call()
  Biplot$Non_Scaled_Data = X
  Biplot$ColumnsMetric=CM
  Biplot$RowsMetric=RM
  
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

  Y= matrixsqrt(CM) %*% X %*% matrixsqrt(RM)
  rownames(Y)=rownames(X)
  colnames(Y)=colnames(X)
  SD = svd(Y, nu = dimension, nv = dimension)
  EV = SD$d^2
  Inertia = round((EV/sum(EV)) * 100, digits = 3)
  CumInertia = cumsum(Inertia)
  a = matrixsqrtinv(CM) %*% SD$u %*% diag(SD$d[1:dimension])
  b = matrixsqrtinv(RM) %*% SD$v %*% diag(SD$d[1:dimension])
  rownames(a) <- RowNames
  colnames(a) <- DimNames
  rownames(b) <- VarNames
  colnames(b) <- DimNames
  CorrXCP = cor(X, a, method = "pearson")
  # Relative contributions of the rows 
  sf = apply((X^2), 1, sum)
  cf=matrix(0,n,dimension)
  for (k in 1:dimension)
    cf[,k]= round((a[,k]^2/ sf)*100, digits = 2)
  rownames(cf) = RowNames
  colnames(cf) = DimNames
  cfacum = t(apply(cf, 1, cumsum))
  # Relative contributions of the columns 
  sc = apply((X^2), 2, sum)
  cc = round(((diag(1/sc)) %*% b^2) * 100, digits = 2)
  rownames(cc) = VarNames
  colnames(cc) = DimNames
  ccacum = t(apply(cc, 1, cumsum))
  if (nfs > 0) {
    as = sup.rows %*% SD$v
    sf = apply((sup.rows^2), 1, sum)
    cfs = round(((diag(1/sf)) %*% as^2) * 100, digits = 2)
    rownames(cfs) = rownames(sup.rows)
    colnames(cfs) = DimNames
  }
  if (ncs > 0) {
    bs = t(sup.cols) %*% SD$u
    sc = apply((sup.cols^2), 2, sum)
    ccs = round(((diag(1/sc)) %*% bs^2) * 100, digits = 2)
    rownames(ccs) = colnames(sup.cols)
    colnames(ccs) = DimNames
  }
  if (alpha <= 1) {
    a = a %*% diag((1/SD$d[1:dimension])^(1 - alpha))
    b = b %*% diag((1/SD$d[1:dimension])^alpha)
    colnames(a) <- DimNames
    colnames(b) <- DimNames
    if (nfs > 0) 
      as = as %*% diag((1/SD$d[1:dimension])^(1 - alpha))
    if (ncs > 0) 
      bs = bs %*% diag((1/SD$d[1:dimension])^alpha)
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
  Biplot$EV = SD$v
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

