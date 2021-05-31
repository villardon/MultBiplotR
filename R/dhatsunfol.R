dhatsunfol <- function(P, D, W = matrix(1, dim(P)[1], dim(P)[2]), modelo, condicion) {
  n = dim(P)[1]
  m = dim(P)[2]
  models = c("Absolute", "Ratio", "Interval", "Ordinal")
  if (condicion == "Matrix") {
    if (modelo == "Identity") {
      DH1 = P
      parametros = matrix(1, m, 1)
    }
    if (modelo == "Ratio") {
      parametros = sum(P * W * D)/sum(P * W * P)
      DH1 = P * parametros
      parametros = parametros * matrix(1, m, 1)
    }
    if (modelo == "Interval") {
      b = sum((P - mean(P)) * W * (D - mean(D)))/sum((P - mean(P)) * W * (P - mean(P)))
      a = mean(D) - b * mean(P)
      DH1 = a + b * P
      parametros = matrix(1, m, 1) %*% c(a, b)
    }
    if (modelo == "Ordinal") {
      rr = isoreg(P, D)
      DH1 = matrix(rr$yf, n, m)
      parametros = matrix(0, m, 1)
    }
  }
  
  if (condicion == "Columns") {
    
    if (modelo == "Absolute") {
      DH1 = P
      parametros = matrix(1, m, 1)
    }
    
    if (modelo == "Ratio") {
      
      parametros = matrix(1, m, 1)
      DH1 = matrix(0, n, m)
      for (j in 1:m) {
        parametros[j] = sum(P[, j] * W[, j] * D[, j])/sum(P[, j] * W[, j] * P[, j])
        parametros[j] = sum(P[, j] * W[, j] * D[, j])/sum(P[, j] * W[, j] * P[, j])
        DH1[, j] = P[, j] * parametros[j]
      }
    }
    if (modelo == "Interval") {
      parametros = matrix(1, m, 2)
      DH1 = matrix(0, n, m)
      for (j in 1:m) {
        parametros[j, 1] = sum((P[, j] - mean(P[, j])) * W[, j] * (D[, j] - mean(D[, j])))/sum((P[, j] - mean(P[, j])) * W[, j] * (P[, j] - mean(P[, 
                                                                                                                                                   j])))
        parametros[j, 2] = mean(D[, j]) - parametros[j, 1] * mean(P[, j])
        DH1[, j] = parametros[j, 1] + P[, j] * parametros[j, 2]
      }
    }
    if (modelo == "Ordinal") {
      DH1 = matrix(0, n, m)
      for (j in 1:m) {
        rr = isoreg(P[, j], D[, j])
        DH1[, j] = rr$yf
        parametros = matrix(0, m, 1)
      }
    }
  }
  Dh = (DH1/sqrt(sum(sum(DH1^2)))) * sqrt(n * m)
  Resul = list(Dh = Dh, Tol = parametros)
  return(Resul)
}

