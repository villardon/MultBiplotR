InitialTransform <- function(X, sup.rows = NULL, sup.cols = NULL, InitTransform="None", transform = "Standardize columns", grouping=NULL) {
  if (is.data.frame(X)) X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  RowNames = rownames(X)
  ColNames = colnames(X)
  
  InitTransforms=c("None", "Log", "Logit")
  if (is.numeric(InitTransform)) 
    InitTransform = InitTransforms[InitTransform]
  
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(transform)) 
    transform = ContinuousDataTransform[transform]

  Data = list()
  Data$InitTransform
  Data$Transformation=transform
  Data$gmean=mean(X, na.rm=TRUE)
  Data$RowMeans=apply(X,1, mean, na.rm=TRUE)
  Data$ColMeans=apply(X,2, mean, na.rm=TRUE)
  Data$RowStdDevs=apply(X,1, sd, na.rm=TRUE)
  Data$ColStdDevs=apply(X,2, sd, na.rm=TRUE)
  Data$ColRanges=apply(X,2, max, na.rm=TRUE)- apply(X,2, min, na.rm=TRUE)
  Data$Grouping=grouping
  
  if (is.null(sup.rows)) 
    nfs = 0
  else {
    nfs = nrow(sup.rows)
    if (!(p == ncol(sup.rows))) 
      stop("The #cols of the supplementary rows must be the same as the #cols of X")
    sup.rows=as.matrix(sup.rows)
  }
  
  if (is.null(sup.cols)) 
    ncs = 0
  else {
    ncs = ncol(sup.cols)
    if (!(n == nrow(sup.cols))) 
      stop("The #rows of the supplementary columns must be the same as the #rows of X")
  }
  
  switch(InitTransform, `Log` = {
    if (sum(which(X<=0)) >0) stop("Initial log transformation is not compatible with negative or zero values")
    X = log(X)
    if (nfs > 0) sup.rows = log(sup.rows)
    if (ncs > 0) sup.cols = log(sup.cols)
  },`Logit` = {
    if (sum(which(X<=0)) >0) stop("Initial logit transformation is not compatible with negative values")
   X= X + 0.01 * (X==0) - 0.01 * (X==1)
   x=log(X/(1-X))
  })
  
  
  nsup = nfs + ncs
  if ((transform == "Double centering") & (nsup > 0)) 
    stop("Double Centering is not compatible with supplementary rows or columns")
  if ((transform == "Divide by the column means and center") & (nsup > 0)) 
    stop("Dividing by the column means and centering is not compatible with supplementary rows or columns")
  if ((transform == "Normalized residuals from independence") & (nsup > 0)) 
    stop("Normalized residuals from independence are not compatible with supplementary rows or columns")
  
  switch(transform, `Substract the global mean` = {
    X = X - Data$gmean
    if (nfs > 0) sup.rows = sup.rows - Data$gmean
    if (ncs > 0) sup.cols = sup.cols - Data$gmean
  }, `Double centering` = {
    X = (diag(n) - matrix(1, n, n)/n) %*% X %*% (diag(p) - matrix(1, p, p)/p)
  }, `Column centering` = {
    X=X- matrix(1,n,1) %*% matrix(Data$ColMeans,1,p)
    if (nfs > 0) sup.rows = (sup.rows - matrix(1, nfs, 1) %*% Data$ColMeans)
    if (ncs > 0) {
      supmeans=apply(sup.cols,2, mean)
      sup.cols=sup.cols- matrix(1,n,1) %*% matrix(supmeans,1,ncs)
      }
  }, `Standardize columns` = {
    X=(X- matrix(1,n,1) %*% matrix(Data$ColMeans,1,p))/(matrix(1,n,1) %*% matrix(Data$ColStdDevs,1,p))
    if (nfs > 0) sup.rows = (sup.rows - matrix(1, nfs, 1) %*% Data$ColMeans)/(matrix(1,nfs,1) %*% matrix(Data$ColStdDevs,1,p))
    if (ncs > 0) {
      supmeans=apply(sup.cols,2, mean)
      supstds=apply(sup.cols,2, sd)
      sup.cols=(sup.cols- matrix(1,n,1) %*% matrix(supmeans,1,ncs))/(matrix(1,n,1) %*% matrix(supstds,1,ncs))
    }
  }, `Row centering` = {
    means = apply(X, 1, mean)
    X = X %*% (diag(p) - matrix(1, p, p)/p)
    if (nfs > 0) {
      means2 = apply(sup.rows, 2, mean)
      for (i in (1:nfs)) {
        sup.rows[i, ] = (sup.rows[i, ] - means2[i])
      }
      if (ncs > 0) {
        for (i in (1:n)) sup.cols[i, ] = (sup.cols[i, ] - sup.cols[i])
      }
    }
  }, `Standardize rows` = {
    means = apply(X, 1, mean)
    stdDevs = apply(X, 1, sd)
    X = solve(diag(stdDevs)) %*% X %*% (diag(p) - matrix(1, p, p)/p)
    if (nfs > 0) {
      means2 = apply(sup.rows, 1, mean)
      stdDEvs2 = apply(sup.rows, 1, sd)
      for (i in (1:nfs)) sup.rows[i, ] = (sup.rows[i, ] - means2[i])/stdDEvs2[i]
    }
    if (ncs > 0) {
      supmeans=apply(sup.cols,2, mean)
      
      for (i in (1:n)) sup.cols[i, ] = (sup.cols[i, ] - means[i])/stdDevs[i]
    }
  }, `Divide by the column means and center` = {
    means = apply(X, 2, mean)
    for (i in (1:p)) X[, i] = X[, i]/ means[i]
    X = (diag(n) - matrix(1, n, n)/n) %*% X
    
  }, `Normalized residuals from independence` = {
    nt = sum(sum(X))
    dr = apply(X,1,sum)
    dc = apply(X,2,sum)
    esp = (t(t(dr)) %*% dc)/nt
    var = t(t(1 - dr/nt)) %*% (1 - dc/nt)
    X = ((X - esp)/sqrt(esp))/sqrt(var)
    
  },`Divide by the range`={
    X=X%*%diag(1/Data$ColRanges)
    
  },`Within groups standardization`={
    if (is.null(grouping)) stop("You need a grouping factor for the within groups standardization")
    if (!is.factor(grouping)) stop("The grouping variable must be a factor")
    g = length(levels(grouping))
    G = Factor2Binary(grouping)
    ng=diag(t(G) %*% G)
    XB= diag(1/ng) %*% t(G) %*% X 
    B = sqrt(diag(ng)) %*% XB
    TSS=apply(X^2,2,sum)
    BSS=apply(B^2,2,sum)
    WSS=(TSS-BSS)/(n-g)
    Data$ColStdDevs=sqrt(WSS)
    X=(X- matrix(1,n,1) %*% matrix(Data$ColMeans,1,p))/(matrix(1,n,1) %*% matrix(WSS,1,p))
  },`Ranks`={
    for (i in 1:p)
      X[,i]=rank(X[,i])
  })
  rownames(X) = RowNames
  colnames(X) = ColNames
  
  Data$X <- X
  Data$sup.rows <- sup.rows
  Data$sup.cols <- sup.cols
  return(Data)
}

