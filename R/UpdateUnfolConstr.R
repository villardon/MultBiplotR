UpdateUnfolConstr <- function(X, Y, W, Dh, ENV, eps=10e-12) {
  n = dim(Dh)[1]
  m = dim(Dh)[2]
  D = DistUnfold(X, Y)
  R = diag(as.vector(W %*% matrix(1, m, 1)))
  C = diag(as.vector(matrix(1, 1, n) %*% W))
  E = matrix(as.numeric(D <= eps), n, m)
  B = W * Dh/(D + E)
  B = B * matrix(as.numeric(E == 0), n, m)
  P = diag(as.vector(B %*% matrix(1, m, 1)))
  Q = diag(as.vector(matrix(1, 1, n) %*% B))
  X1 = P %*% X - B %*% Y
  Y1 = Q %*% Y - t(B) %*% X
  X2 = solve(R) %*% (X1 + W %*% Y)
  X = ENV %*% solve(t(ENV) %*% ENV) %*% t(ENV) %*% X2
  Y = solve(C) %*% (Y1 + t(W) %*% X)
  Conf = list(X = X, Y = Y)
  return(Conf)
}
