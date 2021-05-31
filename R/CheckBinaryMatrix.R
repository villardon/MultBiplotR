# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

CheckBinaryMatrix <- function(x){ 
  for (i in 1:nrow(x))
    for (j in 1:ncol(x))
     if (!((x[i,j]==0) | (x[i,j]==1))) return(FALSE)
  return(TRUE)
}


CheckBinaryVector <- function(x){ 

  for (i in 1:length(x))
      if (!((x[i]==0) | (x[i]==1))) return(FALSE)
  return(TRUE)
}