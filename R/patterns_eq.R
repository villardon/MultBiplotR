patterns_eq <- function(nnodos, dims) {
  I = matrix(1:nnodos)
  for (i in 2:dims) {
    nf = dim(I)[1]
    nc = dim(I)[2]
    I2 = cbind(kronecker(matrix(I[1, ], 1, nc), matrix(1, nnodos, 1)), (1:nnodos))
    for (j in 2:nf) I2 = rbind(I2, cbind(kronecker(matrix(I[j, ], 1, nc), matrix(1, nnodos, 1)), (1:nnodos)))
    I = I2
  }
  return(I)
}