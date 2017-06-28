# MGC stans for Multivariate Gaussian Clustering
MGC <- function(x, NG = 2, init = "km", RemoveOutliers=FALSE, ConfidOutliers=0.995, tolerance = 1e-07, maxiter = 100, show=TRUE, ...) {
  if (is.data.frame(x)) x= as.matrix(x)
  n = dim(x)[1]
  p = dim(x)[2]
  Result=list()
  Result$Data=x
  Result$NG=NG
  # Initial LogLikelihood
  means=matrix(apply(x,2,mean), 1, p)
  xx = x - matrix(1, n, 1) %*% means
  S = (t(xx) %*% xx)/n
  Pr = multnormdens(x, mean = means[1,], sigma = S)
  oldloglik=sum(log(Pr))
  
  Result$InitLoglik=oldloglik
  Result$InitBIC=  2* oldloglik - (p + (p*(p+1)/2)) * log(n)
  
  index=matrix(1,n,NG)
  for (i in 2:NG)
    index[,i]=i
  
  if (init=="km") {
    km = kmeans(x, centers = NG, ...)
    centers = km$centers
    clusters = km$cluster
  }
  
  if (init=="rd"){
    minval=apply(x,2,min)
    maxval=apply(x,2,max)
    rango= maxval - minval
    
    centers = (matrix(runif(NG*p), NG, p) * (matrix(1,NG,1) %*% rango)) + (matrix(1,NG,1) %*% minval)
    D = matrix(0, n, NG)
    for (i in 1:n) for (j in 1:NG) D[i, j] = sqrt(sum((x[i, ] - centers[j, ]) * (x[i, ] - centers[j, ])))
    minim=apply(D,1,min)
    clusters = apply(((D==minim)*index),1,sum)
  }
  # Calculation of the initial variances-covariances matrices 
  C = list()
  alpha = matrix(0, NG, 1)
  for (i in 1:NG) {
    ni = sum(clusters == i)
    xx = x[clusters == i, ] - matrix(1, ni, 1) %*% centers[i, ]
    C[[i]] = (t(xx) %*% xx)/ni
    alpha[i] = sum(clusters == i)/n
  }
  # A posteriori Probabilities for each individual and each cluster
  P = matrix(0, n, NG)
  error = 1
  iter = 0
  
  while ((error > tolerance) & (iter < maxiter)) {
    iter = iter + 1
    for (i in 1:NG) {
      P[, i] = alpha[i] * multnormdens(x, mean = centers[i, ], sigma = C[[i]])
    }
    
    dens=apply(P,1,sum)
    dens=dens/max(dens)
    loglik=sum(log(apply(P,1,sum)))
    error=abs((oldloglik-loglik)/oldloglik)
    oldloglik=loglik
    # Note that dmvnorm could be replaced by any other distribution
    tot = apply(P, 1, sum)
    for (i in 1:NG) {
      P[, i] = P[, i]/tot
    }
    
    OUT=rep(TRUE,n)
    
    if (RemoveOutliers){
      cutoff <- quantile(dens, 1-ConfidOutliers)
      OUT=OUT*(dens>cutoff)
    }
    
    # Updating centers, covariances and pro
    for (i in 1:NG) {
      xx = x * ((P[, i]*OUT) %*% matrix(1, 1, p))
      centers[i, ] = apply(xx, 2, sum)/sum(P[, i]*OUT)
      xx = (x - matrix(1, n, 1) %*% centers[i, ]) * sqrt((P[, i]*OUT) %*% matrix(1, 1, p))
      C[[i]] = (t(xx) %*% xx)/sum(P[, i]*OUT)
    }
    
    # Updating Classification
    MA=apply(P,1,max)
    RE=P>0
    for (i in 1:NG) RE[,i]=(P[,i]==MA)
    Classification=apply(RE*index,1,sum)
    for (i in 1:NG) {
      alpha[i] = sum(Classification == i)/n
    }
    if (show) print(paste(round(iter, digits=0), round(error, digits=9), round(loglik, digits=9)))
  }
  
  
  Result$Centers=centers
  Result$Covariances=C
  Result$P=P
  Result$GroupProbabilities=alpha
  MA=apply(P,1,max)
  RE=P>0
  for (i in 1:NG) RE[,i]=(P[,i]==MA)
  Result$Classification=apply(RE*index,1,sum)
  Result$Outliers=1-OUT
  Result$Loglik=loglik
  Result$BIC = 2*loglik - ((NG-1) + NG* p + (NG * p*(p+1)/2)) * log(n)
  class(Result) <- "MGC"
  return(Result)
}



multnormdens <- function (x, mean = rep(0, p), sigma = diag(p), log = FALSE){
  if (is.vector(x)) 
    x <- matrix(x, ncol = 1)
  p <- ncol(x)
  
  dec <- tryCatch(chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * 
                                                        pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log) 
    logretval
  else exp(logretval)
}


summary.MGC <- function(object, Centers=TRUE, Covariances=TRUE, ...){
  print("MODEL BASED CLUSTERING USING GAUSSIAN MIXTURES")
  print(" No constraints on the covariance Matrices")
  print("___________________________________________________")
  if (Centers) {
    print("Centers")
    print(object$Centers)}
  print("___________________________________________________")
  if (Covariances) {
    print("Covariances")
    print(object$Covariances)}
  print("___________________________________________________")
  print("Sizes of the clusters")
  print(table(object$Classification))
  print("___________________________________________________")
  print(paste("Log-Likelihood : ", object$Loglik))
  print(paste("BIC : ", object$BIC))
}

print.MGC <- function(x, ...){
  print("MODEL BASED CLUSTERING USING GAUSSIAN MIXTURES")
  print(" No constraints on the covariance Matrices")
  print(paste("Log-Likelihood : ", x$Loglik))
  print(paste("BIC : ", x$BIC))
}



