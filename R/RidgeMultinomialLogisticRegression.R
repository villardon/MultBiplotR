# RidgeMultinomialLogisticRegression uses the formula interface but in order to keep
# compatibility with my old programs tha argumant "formula" can also be the response
# and the argument "data" can be the 


RidgeMultinomialLogisticRegression <- function (formula, data, penalization = 0.2, cte = TRUE, tol = 1e-04, 
          maxiter = 200, showIter = FALSE) {
  cl <- match.call()
  if (class(formula)=="formula"){
    mf=model.frame(formula, data=data)
    y=model.extract(mf,"response")
    x=model.matrix(formula, data=data)
    cte=FALSE
  }
  else
  {
    y=formula
    if (!is.factor(y)) stop("The response variable must be a factor")
    x=data
    if (is.matrix(x) & cte) x=cbind(rep(1,n),x)
    if (is.data.frame(x)) x=DataFrame2Matrix4Regression(x, Intercept=cte)
    if (is.vector(x)){
      x=matrix(x, length(x),1)
      if (cte) x=cbind(rep(1,n),x)
    }
  }
  
  catnames=levels(y)
  nlevels=length(catnames)
  Observed=y
  ta=as.matrix(table(y))
  porc=(ta/sum(ta))*100
  itab=cbind(ta, porc)
  colnames(itab)=c("Frequency", "Percentage")
  y=as.integer(y)
  n=nrow(x)
  Model = RidgeMultinomialLogisticFit(y, x, penalization = penalization, tol = tol, 
                     maxiter = maxiter, show = showIter)
  Null = RidgeMultinomialLogisticFit(y, matrix(1, n, 1), penalization = penalization, 
                    tol = tol, maxiter = maxiter, show = showIter)
  Fit = Model
  Fit$call = cl
  rownames(Fit$fitted)=rownames(x)
  colnames(Fit$fitted)=catnames
  Fit$itab=itab
  Fit$Penalization=penalization
  Fit$X=x
  rownames(Fit$Y)=rownames(x)
  colnames(Fit$Y)=catnames
  rownames(Fit$beta)=catnames[-1]
  colnames(Fit$beta)=colnames(x)
  rownames(Fit$stderr)=catnames[-1]
  colnames(Fit$stderr)=colnames(x)
  Fit$NullDeviance = Null$Deviance
  Dif = (Null$Deviance - Model$Deviance)
  Fit$Difference = Dif
  Fit$df = length(Model$beta) - length(Null$beta)
  Fit$p = 1 - pchisq(Dif, df = Fit$df)
  Fit$CoxSnell = 1 - exp(-1 * Dif/n)
  Fit$Nagelkerke = Fit$CoxSnell/(1 - exp((Null$Deviance/(-2)))^(2/n))
  Fit$MacFaden = 1 - (Model$Deviance/Null$Deviance)
  Predicted = rep(0, n)
  correct = rep(0, n)
  for (i in 1:n) {
    Predicted[i]=which(Model$fitted[i, ] == max(Model$fitted[i, ]))
    correct[i] = sum(which(Model$fitted[i, ] == max(Model$fitted[i,])) == y[i])
  }
  if (length(unique(Predicted))==length(catnames)){
    Predicted=factor(Predicted)
  levels(Predicted)=catnames}
  else{
    Predicted=c(Predicted,1:nlevels)
    Predicted=factor(Predicted)
    levels(Predicted)=catnames
    Predicted=Predicted[1:n]
  }
    
  Fit$Table=table(Observed,Predicted)
  Fit$PercentCorrect = sum(correct)/n * 100
  class(Fit)="rmlr"
  return(Fit)
}

RidgeMultinomialLogisticFit <- function (y, x, penalization = 0.2, tol = 1e-04, 
          maxiter = 200, show = FALSE) 
{
  n=length(y)
  if (min(y) == 1) 
    y = y - 1
  p <- ncol(x)
  J = max(y)
  Y = matrix(0, n, J)
  for (i in 1:n) {
    if (y[i] > 0) {
      Y[i, y[i]] = 1
    }
  }
  beta <- matrix(0, J, p)
  yy = as.vector(Y)
  err = 0.1
  iter = 0
  m = matrix(0, n * J, n * J)
  prob = matrix(0, n, J)
  p0 = matrix(0, n, 1)
  xx = kronecker(diag(1, J), x)
  while ((err > tol) & (iter < maxiter)) {
    iter = iter + 1
    betaold = beta
    b = as.vector(t(beta))
    for (i in 1:n) {
      suma = 1
      for (j in 1:J) suma = suma + exp(sum(beta[j, ] * 
                                             x[i, ]))
      for (j in 1:J) prob[i, j] = exp(sum(beta[j, ] * x[i, 
                                                        ]))/suma
      p0[i, 1] = 1/suma
    }
    P = as.vector(prob)
    for (j in 1:J) for (k in 1:J) {
      mj = matrix(0, n, n)
      for (i in 1:n) {
        if (k == j) {
          mj[i, i] = prob[i, j] * (1 - prob[i, j])
        }
        if (!(k == j)) {
          mj[i, i] = -1 * prob[i, j] * prob[i, k]
        }
      }
      m[((j - 1) * n + 1):(j * n), ((k - 1) * n + 1):(k * 
                                                        n)] = mj
    }
    a = (t(xx) %*% m %*% xx) + 2 * diag(J * p) * penalization
    u = (t(xx) %*% matrix((yy - P), n * J, 1)) - 2 * b * 
      penalization
    b = b + solve(a) %*% u
    beta = t(matrix(b, p, J))
    err = sum(sum((betaold - beta) * (betaold - beta)))
  }
  y0 = matrix(0, n, 1)
  for (i in 1:n) if (sum(Y[i, ]) == 0) 
    y0[i] = 1
  cov = solve(a)
  model <- list()
  model$fitted = cbind(p0, prob)
  model$cov = cov
  model$Y = cbind(y0, Y)
  model$beta = beta
  model$stderr = sqrt(t(matrix(diag(cov), p, J)))
  model$logLik = sum(log(model$fitted^model$Y))
  model$Deviance = -2 * sum(log(model$fitted^model$Y))
  model$AIC = model$Deviance + 2 * nrow(beta) * ncol(beta)
  model$BIC = model$Deviance + nrow(beta) * ncol(beta) * log(n)
  class(model) = "rmlr"
  return(model)
}

anova.rmlr <- function(mod1, mod2){
  
}