RidgeBinaryLogisticFit <- function(y, xd, freq, tolerance = 1e-05, maxiter = 100, penalization = 0.2) {
  # Fits a Binary Logistic Regression
  n <- dim(xd)[1]
  m <- dim(xd)[2]
  
  beta = matrix(0, m, 1)
  #freq=freq/sum(freq)
  err = 0.1
  iter = 0
  while ((err > tolerance) & (iter < maxiter)) {
    iter = iter + 1
    betaold = beta
    eta = xd %*% beta
    mu = exp(eta)/(1 + exp(eta))
    v = (mu * (1 - mu)) * freq
    vv = diag(1, n, n)
    diag(vv) <- v
    Imod=diag(m)
    Imod[1,1]=0
    In = (t(xd) %*% vv %*% xd) + 2 * penalization * Imod
    U = t(xd) %*% ((y - mu) * freq) - 2 * penalization * Imod %*% betaold
    beta = betaold + ginv(In) %*% U

    err = sum(abs(betaold - beta))
  }
  return(beta)
}
