indexC1 <- function(P, D) {
  # INDEXC1 Proportion of recovered rank orders 
  #   INDEXC1(P,D) calculates the proportion of rank orders in the vector P 
  #   that are recovered by D.
  
  n = length(P)
  up = 0
  
  for (i in 1:(n - 1)) for (j in (i + 1):n) {
    e = D[i] - D[j]
    work = P[i] - P[j]
    u = 0
    if (work > 0) 
      u = 1
    if (work < 0) 
      u = -1
    work = e * u
    v = 1
    if (work < 0) 
      v = 0
    up = up + v
  }
  
  low = n * (n - 1)/2
  r = up/low
  return(r)
}