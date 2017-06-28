# Update Y values for the genefold algorithm.

UpdateY <- function(SOL, eps = 1e-12) {
  a=SOL$a
  GAMMA=SOL$GAMMA
  X=SOL$X
  Y=SOL$Y
  D=SOL$D
  m = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  vm = matrix(1, 1, m)
  vn = matrix(1, 1, n)
  v = matrix(1, m, m)
  Div = D^(-1)
  Div[which(D < eps)] = 1/eps
  part1 = t(t(vn %*% (((a^2) %*% matrix(1, 1, p)) * X)) %*% vm)
  A = a %*% vm
  wsumGAMMA = (A/m) * (GAMMA %*% v)
  sD = D %*% v
  part2 = matrix(0, m, p)
  part3 = matrix(0, m, p)
  part4 = matrix(0, m, p)
  for (i in 1:p) {
    Xp = X[, i] %*% vm
    Yp = t(vn) %*% t(Y[, i])
    part2i = ((A^2/m) * sD) * ((Yp - Xp) * Div)
    part2[, i] = vn %*% part2i
    part3i = A * (Yp - Xp) * GAMMA * Div
    part3[, i] = vn %*% part3i
    part4i = (wsumGAMMA) * (Xp * Div)
    part4[, i] = vn %*% part4i
  }
  
  Bpart = vn %*% ((wsumGAMMA) * Div)
  A = -2 * (part1 + part2 + part3 + part4)
  B1 = (sum(a^2) + Bpart)
  lambda1 = max(B1)
  int = A + diagonal(2 * B1) %*% Y - 2 * lambda1 * Y
  F2 = -0.5 * int * (lambda1)^(-1)
  Yup = F2
  
  Yup = 2 * Yup - Y
  SOL$Y=Yup
  return(SOL)
}