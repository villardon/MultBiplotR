#   UPDATEGAMMA Monotone regression under a length constraint
#   UpdateGAMMA(DATA,D) performs monotone regression under the constraint  
#   that the variance of GAMMA equals one. It uses the up- and down-
#   block algorithm of Carroll (implemented in the lsqisotonic function of 
#   the stats toolbox), unless the lsqisotonic function results in all equal 
#   values: then a non-negative least squares approach is used.
#   USES the stats toolbox
#   Translated from: K. Van Deun, Department of Psychology, Catholic University of
#   Leuven (BELGIUM)
#   DATA is just a row vector .....

UpdateGAMMA <- function(DATA, D) {
  m = length(DATA)
  n = 1
  iter = 0
  conv = 0
  v = matrix(1, m, m)
  r = lsqisotonic(DATA, D)
  rc = r - r %*% v/m
  
  var = sum(rc^2)
  if (abs(var) < 1e-14) {
    D = t(D)
    J = diag(m) - m^(-1) * v
    index = order(DATA)
    sorted = sort(DATA)
    preforder = order(index)
    x = sort(index)
    G = v
    G[upper.tri(G, diag = FALSE)] = 0
    
    C = matrix(0, m, m)
    for (i in 1:m) C[i, preforder[i]] = 1
    
    S = J %*% C %*% G
    G = t(S) %*% S
    gi = diag(G)
    h = t(S) %*% D
    k = t(D) %*% D + n
    bold = matrix(0, m, 1)
    bold[2] = sqrt(n/gi[2])
    f1 = k - 2 * t(bold) %*% h
    for (i in 3:m) {
      bnew = matrix(0, m, 1)
      bnew[i] = sqrt(n/gi[i])
      f2 = k - 2 * t(bnew) %*% h
      if (f2 - f1 < 0) {
        bold = bnew
        f1 = f2
      }
    }
    bnew = bold
    r = S %*% bnew
    r = t(r)
  } else {
    r = rc/sqrt(var)
  }
  r = r - min(r)
  return(r)
}
