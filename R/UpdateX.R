# Update X values for the genefold algorithm.

UpdateX <- function(SOL, length = 1, lengthconstr = T, eps = 1e-12) {
  # Genefold 
  a=SOL$a
  GAMMA=SOL$GAMMA
  X=SOL$X
  Y=SOL$Y
  D=SOL$D
  m = dim(Y)[1]
  n = dim(X)[1]
  p = dim(X)[2]
  vm = matrix(1, m, 1)
  vn = matrix(1, n, 1)
  v = matrix(1, m, m)
  Div = D^(-1)
  Div[which(D < eps)] = 1/eps
  C3 = Div %*% vm
  A = a %*% t(vm)
  sD = D %*% v
  part1 = matrix(0, n, p)
  part2 = matrix(0, n, p)
  part3 = matrix(0, n, p)
  for (i in 1:p) {
    Xp = X[, i] %*% t(vm)
    Yp = (vn) %*% matrix(Y[, i], 1, m)
    part1i = ((A^2/m) * sD) * ((Xp - Yp) * Div)
    part1[, i] = part1i %*% vm
    part2i = A * (Xp - Yp) * GAMMA * Div
    part2[, i] = part2i %*% vm
    part3i = ((A/m) * (GAMMA %*% v)) * (Yp * Div)
    part3[, i] = part3i %*% vm
  }
  B3 = (a/m) * (GAMMA %*% vm)
  A = -2 * (part1 + part2 + part3)
  BC3 = ((m * (a^2)) + B3 * C3) %*% matrix(1, 1, p)
  int = A
  F2 = -0.5 * int/BC3
  Xup = F2
  Xup = 2 * Xup - X
  if (lengthconstr) {
    sqlength = sqrt(length)
    normX = (Xup * Xup) %*% matrix(1, p, 1)
    minF = -int[which(normX > length), ]
    tr = ((minF * minF) %*% matrix(1, p, 1))^(0.5)
    Xup[which(normX > length), ] = ((sqlength/tr) %*% matrix(1, 1, p)) * minF
  }
  k = which(BC3[, 1] < eps)
  Xup[k, ] = X[k, ]
  SOL$X=Xup
  return(SOL)
}