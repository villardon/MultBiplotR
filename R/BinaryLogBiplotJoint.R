BinaryLogBiplotJoint <- function(x, freq = matrix(1, nrow(x), 1),  dim = 2, ainit=NULL, tolerance = 1e-04, maxiter = 30, penalization=0.2, maxcond = 7, RotVarimax = FALSE, lambda=0.1, ...) {
  # joint algorithm for logistic biplots
  x=as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  
  print("Calculating Initial Row Coordinates")
  
  if (is.null(ainit)) a=svd(x%*%t(x))$u[,1:dim]
  else a=ainit
  
  b = matrix(0, p, dim + 1)
  
  aa=cbind(matrix(1, n, 1), a)
  
  print("Calculating Initial Column Coordinates")
  for (j in 1:p) {
    b[j, ] = RidgeBinaryLogisticFit(x[, j], aa, freq, tolerance, maxiter, penalization)
  }
  
  esp = aa %*% t(b)
  pred = exp(esp)/(1 + exp(esp))
  pred = pred - 1e-04 * as.numeric(pred > 0.9999) + 1e-04 * as.numeric(pred < 1e-04)
  loglikelyhood = sum(sum(x * log(pred) + (1 - x) * log(1 - pred)))
  
  de1 = matrix(0, dim, 1)
  de2 = matrix(0, dim, dim)
  iter = 0
  error = 1
  
  while ((error > tolerance) & (iter < maxiter)) {
    iter = iter + 1
    aold = a
    bold = b
    loglikelyhoodold = loglikelyhood
    
    # Calculating the coordinates of the individuals
    for (i in 1:n) {
      error2 = 1
      iter2 = 0
      separation = 0
      while ((error2 > tolerance) & (iter2 < maxiter) & (separation == 0)) {
        iter2 = iter2 + 1
        a0 = a[i, ]
        mi = b %*% (matrix(c(1, a[i, ])))
        pi = exp(mi)/(1 + exp(mi))
        ww = pi * (1 - pi)
        for (k in 1:dim) de1[k, 1] = -1 * sum(b[, k + 1] * (x[i, ] - pi))
        de1=de1-2*penalization*a0
        for (k in 1:dim) for (l in k:dim) {
          de2[k, l] = -1 * sum(b[, k + 1] * b[, l + 1] * ww)
          de2[l, k] = de2[k, l]
        }
        de2=de2+2*penalization*diag(dim)
        condit=kappa(de2)
        if (condit<maxcond)
          a[i, ] = a0 + (solve(de2) %*% de1)
        error2 = sum(abs(a[i, ] - a0))
        # print(c(iter2, condit,  error2))
      }
    }
    
    aa=cbind(matrix(1, n, 1), a)
    for (j in 1:p) {
      b[j, ] = RidgeBinaryLogisticFit(x[, j], aa, freq, tolerance, maxiter)
    }
    
    esp = aa %*% t(b)
    pred = exp(esp)/(1 + exp(esp))
    pred = pred - 1e-04 * as.numeric(pred > 0.9999) + 1e-04 * as.numeric(pred < 1e-04)
    loglikelyhood = sum(sum(x * log(pred) + (1 - x) * log(1 - pred)))
    
    error = (loglikelyhood -loglikelyhoodold )/abs(loglikelyhood)
    if (error < 0){
      a=aold
      b=bold}
    else {
      aold = a
    bold = b}
    print(c(iter, loglikelyhood, error))
  }
  
  if (RotVarimax) {
    varimax(b)
    B = varimax(b[, 1:dim + 1], normalize = FALSE)
    a = a %*% B$rotmat
    b[, 1:dim + 1] = b[, 1:dim + 1] %*% B$rotmat
    esp = cbind(freq, a) %*% t(b)
    pred = exp(esp)/(1 + exp(esp))
  }
  
  pred2 = matrix(as.numeric(pred > 0.5), n, p)
  acier = matrix(as.numeric(round(x) == pred2), n, p)
  acierfil = apply(acier,1,sum)
  aciercol = apply(acier,2,sum)
  
  gfit = (sum(sum(acier))/(n * p)) * 100
  
  
  pred = pred - 1e-04 * as.numeric(pred > 0.9999) + 1e-04 * as.numeric(pred < 1e-04)
  
  esp0 = freq %*% b[, 1]
  pred0 = exp(esp0)/(1 + exp(esp0))
  
  d1 = -2 * apply(x * log(pred0) + (1 - x) * log(1 - pred0),2,sum)
  d2 = -2 * apply(x * log(pred) + (1 - x) * log(1 - pred),2,sum)
  
  d = d1 - d2
  ps = matrix(0, p, 1)
  for (j in 1:p) ps[j] = 1 - pchisq(d[j], 1)
  
  resid = x - pred
  r2 = 1 - (sum(sum(resid^2))/sum(sum(x^2)))
  pred2 = (esp > 0.5)
  corrfil = (acierfil/p) * 100
  corrcol = (aciercol/n) * 100
  rownames(b)=colnames(x)
  result=list(RowCoordinates=a, ColumnParameters=b, Deviances=d, pValues=ps, GlobalPercent= gfit, RowPercents=corrfil, ColPercents=corrcol, R2=r2)
  class(result) = "Binary.Logistic.Biplot"
  return(result)
}

