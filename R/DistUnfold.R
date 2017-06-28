DistUnfold <- function(X, Y) {
  n = dim(X)[1]
  m = dim(Y)[1]
  D = matrix(0, n, m)
  for (i in 1:n) for (j in 1:m) D[i, j] = sqrt(sum((X[i, ] - Y[j, ]) * (X[i, ] - Y[j, ])))
  return(D)
}