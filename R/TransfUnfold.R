TransfUnfold <- function(f, trans, offset) {
  n = dim(f)[1]
  m = dim(f)[2]
  transformations=c("None", "Gaussian", "Column Percent", "Gaussian Columns", "Inverse Square Root", "Divide by Column Maximum")
  if (is.numeric(trans)) trans=transformations[trans]
  if (trans == "None") 
    P = f
  if (trans == "Gaussian") {
    f= f + matrix(as.numeric(f == 0), n, m) * offset
    maxco = max(f)
    P = sqrt(-2 * log((f + matrix(as.numeric(f == 0), n, m) * offset)/maxco))
  }
  if (trans == "Column Percent") {
    sumco = apply(f, 2, sum)
    P = f %*% solve(diag(sumco))
  }
  if (trans == "Gaussian Columns") {
    f= as.matrix(f + matrix(as.numeric(f == 0), n, m) * offset)
    maxco = apply(f, 2, max)
    P =sqrt(-2 * log(f %*% diag(1/maxco)))
  }
  if (trans == "Inverse Square Root") {
    P = 1/sqrt((f + matrix(as.numeric(f == 0), n, m) * offset))
  }
  if (trans == "Divide by Column Maximum") {
    maxco = apply(f, 2, max)
    P = f %*% solve(diag(maxco))
  }
  colnames(P) = colnames(f)
  rownames(P) = rownames(f)		
  return(P)
}