matrixsqrtinv <- function(S, tol = sqrt(.Machine$double.eps))
{
  ## Square root of a Matrix
  s <- svd(S)
  nz <- s$d > tol
  S12=s$u[, nz] %*% diag(1/sqrt(s$d[nz])) %*% t(s$v[, nz])
  return(S12)
}