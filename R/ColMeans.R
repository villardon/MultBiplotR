ColMeans <- function(X) {
  dimens = dim(X)
  n = dimens[1]
  p = dimens[2]
  Means = rep(0, p)
  for (i in (1:p)) Means[i] = mean(X[, i])
  return(Means)
}