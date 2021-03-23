# Using non vectorial programming (for rather than matrices) is less efficient
# except for large data matices
TransformIni <- function(X, InitTransform="None", transform = "Standardize columns") {
  n = nrow(X)
  p = ncol(X)
  RowNames = rownames(X)
  ColNames = colnames(X)
  
  InitTransforms=c("None", "Log", "Logit")
  if (is.numeric(InitTransform)) 
    InitTransform = InitTransforms[InitTransform]
  
  
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center", "Normalized residuals from independence")
  if (is.numeric(transform)) 
    transform = ContinuousDataTransform[transform]
  
  switch(InitTransform, `Log` = {
    if (sum(which(X<=0)) >0) stop("Initial log transformation is not compatible with negative or zero values")
    X = log(X)
  },`Logit` = {
    if (sum(which(X<=0)) >0) stop("Initial logit transformation is not compatible with negative values")
    X= X + 0.01 * (X==0) - 0.01 * (X==1)
    x=log(X/(1-X))
  })
  
  
  switch(transform, `Substract the global mean` = {
    gmean = mean(X)
    X = X - gmean
  }, `Double centering` = {
    X = (diag(n) - matrix(1, n, n)/n) %*% X %*% (diag(p) - matrix(1, p, p)/p)
  }, `Column centering` = {
    means = apply(X, 2, mean)
    X=X- matrix(1,n,1) %*% matrix(means,1,p)
  }, `Standardize columns` = {
    means = apply(X, 2, mean)
    stdDevs = apply(X, 2, sd)
    X=(X- matrix(1,n,1) %*% matrix(means,1,p))/(matrix(1,n,1) %*% matrix(stdDevs,1,p))
  }, `Row centering` = {
    means = apply(X, 1, mean)
    X = X %*% (diag(p) - matrix(1, p, p)/p)
  }, `Standardize rows` = {
    means = apply(X, 1, mean)
    stdDevs = apply(X, 1, sd)
    X = solve(diag(stdDevs)) %*% X %*% (diag(p) - matrix(1, p, p)/p)
  }, `Divide by the column means and center` = {
    means = apply(X, 2, mean)
    for (i in (1:p)) X[, i] = X[, i]/means[i]
    X = (diag(n) - matrix(1, n, n)/n) %*% X
  }, `Normalized residuals from independence` = {
    nt = sum(sum(X))
    dr = apply(X,1,sum)
    dc = apply(X,2,sum)
    esp = (t(t(dr)) %*% dc)/nt
    var = t(t(1 - dr/nt)) %*% (1 - dc/nt)
    X = ((X - esp)/sqrt(esp))/sqrt(var)
  },`Divide by the range`={
    Rangos=apply(X,2,max)-apply(X,2,min)
    X=X%*%diag(1/Rangos)
  })
  rownames(X) = RowNames
  colnames(X) = ColNames
  return(X)
}
